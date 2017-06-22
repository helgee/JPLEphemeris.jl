using AstroDynBase
using LightGraphs

import AstroDynBase: position, velocity, state, position!, velocity!, state!

export SPK, position, velocity, state, position!, velocity!, state!,
    segments, print_segments

const SECONDS_PER_DAY = 86400
const SIZE_FLOAT64 = sizeof(Float64)

type OutOfRangeError <: Exception
    date::Float64
    startdate::Float64
    finaldate::Float64
end

Base.showerror(io::IO, err::OutOfRangeError) = print(io, "The requested date $(err.date) is outside the intervall ($(err.startdate), $(err.finaldate)).")

type Segment
    name::String
    firstsec::Float64
    lastsec::Float64
    firstdate::Float64
    lastdate::Float64
    target::Int
    center::Int
    frame::Int
    spktype::Int
    firstaddr::Int
    lastaddr::Int
    firstword::Int
    lastword::Int
    initialsecond::Float64
    intlen::Float64
    rsize::Int
    n_records::Int
    order::Int
    cached_record::Int
    cache::Matrix{Float64}
end

jd(sec) = 2451545 + sec/SECONDS_PER_DAY
seconds(jd) = (jd - 2451545)*SECONDS_PER_DAY

function Segment(daf, name, record)
    firstsec, lastsec = reinterpret_getindex(Float64, record, (1, 9), daf.little)
    target, center, frame, spktype, firstaddr, lastaddr =
        reinterpret_getindex(Int32, record, (17, 21, 25, 29, 33, 37), daf.little)
    if spktype != 2
        error("Type $spktype SPK file detected. Only Type 2 SPK files are supported.")
    end
    i0 = lastaddr * SIZE_FLOAT64 - 4 * SIZE_FLOAT64 + 1
    init, intlen, rsize, n_records =
        reinterpret_getindex(Float64, daf.array, (i0, i0 + 8, i0 + 16, i0 + 24), daf.little)
    n_records = round(Int32, n_records)
    order = Int((rsize - 2) ÷ 3)
    Segment(
        name,
        firstsec,
        lastsec,
        jd(firstsec),
        jd(lastsec),
        target,
        center,
        frame,
        spktype,
        firstaddr,
        lastaddr,
        firstaddr*SIZE_FLOAT64 - SIZE_FLOAT64 + 1,
        lastaddr*SIZE_FLOAT64 - SIZE_FLOAT64*4,
        init,
        intlen,
        round(Int32, rsize),
        n_records,
        order,
        -1,
        zeros(3, order),
    )
end

type SPK <: Ephemeris
    daf::DAF
    segments::Dict{Int,Dict{Int,Segment}}
    paths::Dict{Int,Dict{Int,Vector{Int}}}
end

Base.show(io::IO, spk::SPK) = print(io, "SPK($(spk.segments[0][1].name))")

function SPK(filename)
    daf = DAF(filename)
    segments = Dict{Int, Dict{Int, Segment}}()
    graph = Graph()
    to_id = Dict{Int,Int}()
    to_body = Dict{Int,Int}()
    for (name, summary) in getsummaries(daf)
        seg = Segment(daf, name, summary)
        if haskey(segments, seg.center)
            merge!(segments[seg.center], Dict(seg.target=>seg))
        else
            merge!(segments, Dict(seg.center=>Dict(seg.target=>seg)))
        end

        if !(seg.center in keys(to_id))
            add_vertex!(graph)
            merge!(to_id, Dict(seg.center=>nv(graph)))
            merge!(to_body, Dict(nv(graph)=>seg.center))
        end
        if !(seg.target in keys(to_id))
            add_vertex!(graph)
            merge!(to_id, Dict(seg.target=>nv(graph)))
            merge!(to_body, Dict(nv(graph)=>seg.target))
        end
        add_edge!(graph, to_id[seg.center], to_id[seg.target])
    end

    paths = Dict{Int,Dict{Int,Vector{Int}}}()
    for (origin, oid) in to_id
        d = dijkstra_shortest_paths(graph, oid)
        for (target, tid) in to_id
            origin == target && continue
            path = map(x->to_body[x], enumerate_paths(d, tid))
            if haskey(paths, origin)
                merge!(paths[origin], Dict(target=>path))
            else
                merge!(paths, Dict(origin=>Dict(target=>path)))
            end
        end
    end

    SPK(daf, segments, paths)
end

segments(spk::SPK) = spk.segments

function list_segments(spk::SPK)
    s = String[]
    for (k,v) in spk.segments
        for l in keys(v)
            push!(s, "$(name_from_naifid(k)) ($k) => $(name_from_naifid(l)) ($l)")
        end
    end
    return sort!(s, lt=segstrlt)
end

function print_segments(spk::SPK)
    s = list_segments(spk)
    println(join(s, "\n"))
end

function segstrlt(a::String, b::String)
   rex = r"\([0-9]*\)$"
   ma = match(rex, a)
   mb = match(rex, b)
   ia = parse(Int, a[ma.offset+1:end-1])
   ib = parse(Int, b[mb.offset+1:end-1])
   return ia < ib
end

@inline function checkdate(seg::Segment, tdb::Float64)
    if !(seg.firstdate <= tdb <= seg.lastdate)
        throw(OutOfRangeError(tdb, seg.firstdate, seg.lastdate))
    end
end

@inline function getrecordnum(seg, tdb, tdb2)
    checkdate(seg, tdb+tdb2)
    secs = (seconds(tdb) - seg.initialsecond) + tdb2 * SECONDS_PER_DAY
    recordnum, frac = divrem(secs, seg.intlen)
    recordnum = round(Int, recordnum)
    if recordnum == seg.n_records
        recordnum -= 1
        frac = seg.intlen
    end
    recordnum, frac
end

@inline function update_cache!(spk::SPK, seg::Segment, recordnum)
    components = 3
    seg.cached_record = recordnum
    # Drop the MID and RADIUS values
    first = seg.firstword + SIZE_FLOAT64 * seg.rsize * recordnum + SIZE_FLOAT64 * 2
    ptr = Ptr{Float64}(pointer(spk.daf.array, first))

    cache = unsafe_wrap(Array, Ptr{Float64}(ptr), (seg.order, components), false)
    if !spk.daf.little
        transpose!(seg.cache, ntoh.(copy(cache)))
    else
        transpose!(seg.cache, cache)
    end
end

@inline function chebyshev(order, intlen, frac)
    x = Array{Float64}(order)
    x[1] = 1.0
    x[2] = 2.0 * frac/intlen - 1.0
    @inbounds for i = 3:order
        x[i] = 2.0 * x[2] * x[i-1] - x[i-2]
    end
    x
end

@inline function chebyshev_deriv(x, order, intlen)
    t = zeros(Float64, order)
    t[2] = 1.0
    if order > 2
        t[3] = 4.0 * x[2]
        @inbounds for i = 4:order
            t[i] = 2.0 * x[2] *t[i-1] - t[i-2] + x[i-1] + x[i-1]
        end
    end
    t .*= 2.0
    t ./= intlen
end

@inline function position!(r, x::AbstractArray, seg::Segment, sign::Float64)
    @inbounds @simd for i = 1:3
        for j = 1:seg.order
            r[i] += sign * seg.cache[i, j] * x[j]
        end
    end
    r
end

@inline function position!(r, spk::SPK, seg::Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    x = chebyshev(seg.order, seg.intlen, frac)
    position!(r, x, seg, sign)
end

@inline function velocity!(v, x::AbstractArray, seg::Segment, sign::Float64)
    t = chebyshev_deriv(x, seg.order, seg.intlen)
    @inbounds @simd for i = 1:3
        for j = 1:seg.order
            v[i] += sign * seg.cache[i, j] * t[j]
        end
    end
    v
end

@inline function velocity!(v, spk::SPK, seg::Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    x = chebyshev(seg.order, seg.intlen, frac)
    velocity!(v, x, seg, sign)
end


@inline function state!(s, spk::SPK, seg::Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    x = chebyshev(seg.order, seg.intlen, frac)
    @views begin
        position!(s[1:3], x, seg, sign)
        velocity!(s[4:6], x, seg, sign)
    end
    s
end

@inline function findsegment(segments, origin, target)
    if !(origin in keys(segments) || target in keys(segments))
        error("No segment '$origin'->'$target' available.")
    end
    sign = 1.0
    if target < origin
        origin, target = target, origin
        sign = -1.0
    end
    return segments[origin][target], sign
end

for (f, n) in zip((:state, :velocity, :position), (6, 3, 3))
    fmut = Symbol(f, "!")
    @eval begin
        function $fmut(arr, spk::SPK, ep::TDBEpoch, from::Type{C1}, to::Type{C2}) where {C1<:CelestialBody, C2<:CelestialBody}
            path = spk.paths[naif_id(from)][naif_id(to)]
            jd1 = julian1(ep)
            jd2 = julian2(ep)

            $fmut(arr, spk, path[1], path[2], jd1, jd2)
            for (origin, target) in zip(path[2:end-1], path[3:end])
                $fmut(arr, spk, origin, target, jd1, jd2)
            end
            arr
        end

        function $f(spk::SPK, ep::TDBEpoch, from::Type{C1}, to::Type{C2}) where {C1<:CelestialBody, C2<:CelestialBody}
            arr = zeros($n)
            $fmut(arr, spk, ep, from, to)
        end

        function $fmut(arr, spk::SPK, center::Int, target::Int, tdb::Float64, tdb2::Float64=0.0)
            seg, sign = findsegment(spk.segments, center, target)
            $fmut(arr, spk, seg, sign, tdb, tdb2)
        end

        function $fmut(arr, spk::SPK, target::Int, tdb::Float64, tdb2::Float64=0.0)
            seg = spk.segments[0][target]
            $fmut(arr, spk, seg, 1.0, tdb, tdb2)
        end

        function $fmut(arr, spk::SPK, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            $fmut(arr, spk, naifid(target), tdb, tdb2)
        end

        function $fmut(arr, spk::SPK, center::AbstractString, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            $fmut(arr, spk, naifid(center), naifid(target), tdb, tdb2)
        end

        function $fmut(arr, spk::SPK, center::Int, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            $fmut(arr, spk, center, naifid(target), tdb, tdb2)
        end

        function $fmut(arr, spk::SPK, center::AbstractString, target::Int, tdb::Float64, tdb2::Float64=0.0)
            $fmut(arr, spk, naifid(center), target, tdb, tdb2)
        end

        $f(spk::SPK, target, tdb::Float64, tdb2::Float64=0.0) =
            $fmut(zeros($n), spk, target, tdb, tdb2)

        $f(spk::SPK, center, target, tdb::Float64, tdb2::Float64=0.0) =
            $fmut(zeros($n), spk, center, target, tdb, tdb2)
    end
end
