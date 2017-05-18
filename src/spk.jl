using AstroDynBase
import AstroDynBase: position, velocity, state

export SPK, position, velocity, state, segments, print_segments

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
    target::Int32
    center::Int32
    frame::Int32
    spktype::Int32
    firstaddr::Int32
    lastaddr::Int32
    firstword::Int32
    lastword::Int32
    initialsecond::Float64
    intlen::Float64
    rsize::Int32
    n_records::Int32
    cache::Matrix{Float64}
    cached_record::Int
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
        Matrix{Float64}(0,0),
        -1,
    )
end

type SPK <: Ephemeris
    daf::DAF
    segments::Dict{Int, Dict{Int, Segment}}
end

Base.show(io::IO, spk::SPK) = print(io, "SPK($(spk.segments[0][1].name))")

function SPK(filename)
    daf = DAF(filename)
    segments = Dict{Int, Dict{Int, Segment}}()
    for (name, summary) in getsummaries(daf)
        seg = Segment(daf, name, summary)
        if haskey(segments, seg.center)
            merge!(segments[seg.center], Dict(seg.target=>seg))
        else
            merge!(segments, Dict(seg.center=>Dict(seg.target=>seg)))
        end
    end
    SPK(daf, segments)
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

function checkdate(seg::Segment, tdb::Float64)
    if !(seg.firstdate <= tdb <= seg.lastdate)
        throw(OutOfRangeError(tdb, seg.firstdate, seg.lastdate))
    end
end

function getcoefficients(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    checkdate(seg, tdb+tdb2)
    components = 3
    order = (seg.rsize - 2) ÷ 3
    secs = (seconds(tdb) - seg.initialsecond) + tdb2*SECONDS_PER_DAY
    recordnum, frac = divrem(secs, seg.intlen)
    recordnum = round(Int, recordnum)
    if recordnum == seg.n_records
        recordnum -= 1
        frac = seg.intlen
    end
    if recordnum == seg.cached_record
        c = seg.cache
    else
        # Drop the MID and RADIUS values
        first = seg.firstword + SIZE_FLOAT64*seg.rsize*recordnum + SIZE_FLOAT64*2
        last = seg.firstword + SIZE_FLOAT64*seg.rsize*(recordnum+1) - 1
        cv = reinterpret_getindex.(Float64, (spk.daf.array,), first:8:last, spk.daf.little)
        c = reshape(cv, (order, components))'
        seg.cache = c
        seg.cached_record = recordnum
    end
    x = Array{Float64}(order)
    tc = 2.0 * frac/seg.intlen - 1.0
    x[1] = 1.0
    x[2] = tc
    twotc = tc + tc
    @inbounds for i = 3:order
        x[i] = twotc*x[i-1] - x[i-2]
    end
    c, x, seg.intlen, twotc, order
end

function position(c::Matrix, x::Vector, order::Int)
    r = zeros(3)
    @inbounds @simd for i = 1:3
        for j = 1:order
            r[i] += c[i, j] * x[j]
        end
    end
    return r
end

function position(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc, order = getcoefficients(spk, seg, tdb, tdb2)
    position(c, x, order)
end

function velocity(c::Matrix, x::Vector, dt::Float64, twotc::Float64, order::Int)
    v = zeros(Float64, 3)
    t = zeros(Float64, order)
    t[2] = 1.0
    if order > 2
        t[3] = twotc + twotc
        for i = 4:order
            t[i] = twotc*t[i-1] - t[i-2] + x[i-1] + x[i-1]
        end
    end
    t *= 2.0
    t /= dt
    @inbounds @simd for i = 1:3
        for j = 1:order
            v[i] += c[i, j] * t[j]
        end
    end
    return v
end

function velocity(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc, order = getcoefficients(spk, seg, tdb, tdb2)
    velocity(c, x, dt, twotc, order)
end


function state(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc, order = getcoefficients(spk, seg, tdb, tdb2)
    r = position(c, x, order)
    v = velocity(c, x, dt, twotc, order)
    r, v
end

function findsegment(segments, origin, target)
    if origin in keys(segments) && target in keys(segments[origin])
        factor = 1.0
        return segments[origin][target], factor
    elseif target in keys(segments) && origin in keys(segments[target])
        factor = -1.0
        return segments[target][origin], factor
    else
        error("No segment '$origin'->'$target' available.")
    end
end

function findpath(origin, target)
    if target == parent(origin) || parent(target) == origin
        return [origin, target]
    elseif target == parent(parent(origin))
        return [origin, parent(origin), target]
    elseif origin == parent(parent(target))
        return [origin, parent(target), target]
    elseif parent(target) == parent(origin)
        return [origin, parent(origin), target]
    elseif parent(parent(target)) == parent(origin) ||
        parent(target) == parent(parent(origin))
        return [origin, parent(origin), parent(target), target]
    elseif parent(parent(target)) == parent(parent(origin))
        return [origin, parent(origin), parent(parent(origin)), parent(target), target]
    end
end

for (f, u) in zip((:state, :velocity, :position), ((:km, :kps), (:kps,), (:km,)))
    @eval begin
        function $f(spk::SPK, ep::TDBEpoch, from::Type{C1}, to::Type{C2}) where {C1<:CelestialBody, C2<:CelestialBody}
            path = findpath(from, to)
            jd1 = julian1_strip(ep)
            jd2 = julian2_strip(ep)
            length(path) == 2 && $f(spk, naif_id(from), naif_id(to), jd1, jd2)

            res = $f(spk, naif_id(path[1]), naif_id(path[2]), jd1, jd2)
            for (origin, target) in zip(path[2:end-1], path[3:end])
                res = res .+ $f(spk, naif_id(origin), naif_id(target), jd1, jd2)
            end
            res .* ($(u...),)
        end

        function $f(spk::SPK, center::Int, target::Int, tdb::Float64, tdb2::Float64=0.0)
            seg, factor = findsegment(spk.segments, center, target)
            $f(spk, seg, tdb, tdb2) .* factor
        end

        function $f(spk::SPK, target::Int, tdb::Float64, tdb2::Float64=0.0)
            seg = spk.segments[0][target]
            $f(spk, seg, tdb, tdb2)
        end

        function $f(spk::SPK, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            $f(spk, naifid(target), tdb, tdb2)
        end

        function $f(spk::SPK, center::AbstractString, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            $f(spk, naifid(center), naifid(target), tdb, tdb2)
        end

        function $f(spk::SPK, center::Int, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            $f(spk, center, naifid(target), tdb, tdb2)
        end

        function $f(spk::SPK, center::AbstractString, target::Int, tdb::Float64, tdb2::Float64=0.0)
            $f(spk, naifid(center), target, tdb, tdb2)
        end
    end
end
