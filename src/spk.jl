using LinearAlgebra: transpose!
using AstroBase: NAIFId, TDBEpoch, from_naifid, julian_twopart, value, SECONDS_PER_DAY

import AstroBase.Interfaces:
    AbstractEphemeris,
    position,
    position!,
    velocity,
    velocity!,
    position_velocity,
    position_velocity!

export SPK,
    segments,
    print_segments,
    position,
    position!,
    velocity,
    velocity!,
    position_velocity,
    position_velocity!

const SIZE_FLOAT64 = sizeof(Float64)

struct OutOfRangeError <: Exception
    date::Float64
    startdate::Float64
    finaldate::Float64
end

Base.showerror(io::IO, err::OutOfRangeError) = print(io,
   "The requested date $(err.date) is outside the interval ($(err.startdate), $(err.finaldate)).")

mutable struct Segment
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
    x::Vector{Float64}
    t::Vector{Float64}
end

jd(sec) = 2451545.0 + sec / SECONDS_PER_DAY
seconds(jd) = (jd - 2451545.0) * SECONDS_PER_DAY

function Segment(daf, name, record)
    firstsec, lastsec = reinterpret_getindex(Float64, record, (1, 9), daf.little)
    target, center, frame, spktype, firstaddr, lastaddr =
        reinterpret_getindex(Int32, record, (17, 21, 25, 29, 33, 37), daf.little)
    if spktype != 2
        throw(ArgumentError("Type $spktype SPK file detected. Only Type 2 SPK files are supported."))
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
        zeros(order),
        zeros(order),
    )
end

struct SPK <: AbstractEphemeris
    daf::DAF
    segments::Dict{Int,Dict{Int,Segment}}
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
            push!(s, "$(from_naifid(k)) ($k) => $(from_naifid(l)) ($l)")
        end
    end
    sort!(s, lt=segstrlt)
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
   ia < ib
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

    cache = unsafe_wrap(Array, Ptr{Float64}(ptr), (seg.order, components), own=false)
    if !spk.daf.little
        transpose!(seg.cache, ntoh.(copy(cache)))
    else
        transpose!(seg.cache, cache)
    end
end

@inline function chebyshev!(seg, frac)
    seg.x[1] = 1.0
    seg.x[2] = 2.0 * frac / seg.intlen - 1.0
    @inbounds for i = 3:seg.order
        seg.x[i] = 2.0 * seg.x[2] * seg.x[i-1] - seg.x[i-2]
    end
end

@inline function chebyshev_deriv!(seg)
    seg.t[2] = 1.0
    if seg.order > 2
        seg.t[3] = 4.0 * seg.x[2]
        @inbounds for i = 4:seg.order
            seg.t[i] = 2.0 * seg.x[2] * seg.t[i-1] - seg.t[i-2] +
                seg.x[i-1] + seg.x[i-1]
        end
    end
    seg.t .*= 2.0
    seg.t ./= seg.intlen
end

@inline function position!(r, seg::Segment, sign::Float64)
    @inbounds @simd for i = 1:3
        for j = 1:seg.order
            r[i] += sign * seg.cache[i, j] * seg.x[j]
        end
    end
    r
end

@inline function position!(r, spk::SPK, seg::Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    chebyshev!(seg, frac)
    position!(r, seg, sign)
end

@inline function velocity!(v, seg::Segment, sign::Float64)
    chebyshev_deriv!(seg)
    @inbounds @simd for i = 1:3
        for j = 1:seg.order
            v[i] += sign * seg.cache[i, j] * seg.t[j]
        end
    end
    v
end

@inline function velocity!(v, spk::SPK, seg::Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    chebyshev!(seg, frac)
    velocity!(v, seg, sign)
end


@inline function position_velocity!(pos,
                                    vel,
                                    spk::SPK,
                                    seg::Segment,
                                    sign::Float64,
                                    tdb::Float64,
                                    tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    chebyshev!(seg, frac)
    position!(pos, seg, sign)
    velocity!(vel, seg, sign)
    pos, vel
end

@inline function findsegment(segments, origin, target)
    if !(origin in keys(segments) || target in keys(segments))
        throw(ArgumentError("No segment '$origin'->'$target' available."))
    end
    sign = 1.0
    if target < origin
        origin, target = target, origin
        sign = -1.0
    end
    segments[origin][target], sign
end


function position!(pos, spk::SPK, ep::TDBEpoch, from::NAIFId, to::NAIFId)
    seg, sign = findsegment(spk.segments, from, to)
    jd1, jd2 = value.(julian_twopart(ep))
    position!(pos, spk, seg, sign, jd1, jd2)
end

function velocity!(vel, spk::SPK, ep::TDBEpoch, from::NAIFId, to::NAIFId)
    seg, sign = findsegment(spk.segments, from, to)
    jd1, jd2 = value.(julian_twopart(ep))
    velocity!(vel, spk, seg, sign, jd1, jd2)
end

function position_velocity!(pos, vel, spk::SPK, ep::TDBEpoch, from::NAIFId, to::NAIFId)
    seg, sign = findsegment(spk.segments, from, to)
    jd1, jd2 = value.(julian_twopart(ep))
    position_velocity!(pos, vel, spk, seg, sign, jd1, jd2)
end

