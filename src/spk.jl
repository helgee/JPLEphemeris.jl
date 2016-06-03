import Base.position

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
    name::ASCIIString
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
    firstsec, lastsec = reinterpret(Float64, record[1:16], daf.little)
    target, center, frame, spktype, firstaddr, lastaddr = reinterpret(Int32, record[17:end], daf.little)
    if spktype != 2
        error("Type $spktype SPK file detected. Only Type 2 SPK files are supported.")
    end
    init, intlen, rsize, n_records = reinterpret(Float64, daf.array[lastaddr*SIZE_FLOAT64-4*SIZE_FLOAT64+1:lastaddr*SIZE_FLOAT64], daf.little)
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
        Matrix{Float64}(),
        -1,
    )
end

type SPK
    daf::DAF
    segments::Dict{Int, Dict{Int, Segment}}
end

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
    s = ASCIIString[]
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

function segstrlt(a::ASCIIString, b::ASCIIString)
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
        last = seg.firstword + SIZE_FLOAT64*seg.rsize*(recordnum+1)
        cv = reinterpret(Float64, spk.daf.array[first:last], spk.daf.little)
        c = reshape(cv, (order, components))'
        seg.cache = c
        seg.cached_record = recordnum
    end
    x = Array(Float64, order)
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
    [r;v]
end

for (f, n) in zip((:state, :velocity, :position), (6, 3, 3))
    @eval begin
        function ($f)(spk::SPK, center::Int, target::Int, tdb::Float64, tdb2::Float64=0.0)
            seg = spk.segments[center][target]
            ($f)(spk, seg, tdb, tdb2)
        end

        function ($f)(spk::SPK, target::Int, tdb::Float64, tdb2::Float64=0.0)
            seg = spk.segments[0][target]
            ($f)(spk, seg, tdb, tdb2)
        end

        function ($f)(spk::SPK, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            ($f)(spk, naifid(target), tdb, tdb2)
        end

        function ($f)(spk::SPK, center::AbstractString, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            ($f)(spk, naifid(center), naifid(target), tdb, tdb2)
        end

        function ($f)(spk::SPK, center::Int, target::AbstractString, tdb::Float64, tdb2::Float64=0.0)
            ($f)(spk, center, naifid(target), tdb, tdb2)
        end

        function ($f)(spk::SPK, center::AbstractString, target::Int, tdb::Float64, tdb2::Float64=0.0)
            ($f)(spk, naifid(center), target, tdb, tdb2)
        end

        function ($f)(spk::SPK, center, target, tdb::AbstractArray, tdb2::AbstractArray=zeros(length(tdb)))
            m = length(tdb)
            if m != length(tdb2)
                error("'tdb' and 'tdb2' must have the same length.")
            end
            out = zeros(m, $n)
            for (i, t, t2) in zip(1:m, tdb, tdb2)
                out[i,:] = ($f)(spk, center, target, t, t2)
            end
            return out
        end

        function ($f)(spk::SPK, target, tdb::AbstractArray, tdb2::AbstractArray=zeros(length(tdb)))
            ($f)(spk, 0, target, tdb, tdb2)
        end
    end
end
