import Base.position

export SPK, position, velocity, state

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
    n::Int32
end

jd(sec) = 2451545 + sec/SECONDS_PER_DAY
seconds(jd) = (jd - 2451545)*SECONDS_PER_DAY

function Segment(daf, name, record)
    firstsec, lastsec = reinterpret(Float64, record[1:16], daf.little)
    target, center, frame, spktype, firstaddr, lastaddr = reinterpret(Int32, record[17:end], daf.little)
    if spktype != 2
        error("Only Type 2 SPK files are supported.")
    end
    init, intlen, rsize, n = reinterpret(Float64, daf.array[lastaddr*SIZE_FLOAT64-4*SIZE_FLOAT64+1:lastaddr*SIZE_FLOAT64], daf.little)
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
        round(Int32, n),
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
    if recordnum == seg.n
        seg -= 1
        frac = seg.intlen
    end
    # Drop the MID and RADIUS values
    first = seg.firstword + SIZE_FLOAT64*seg.rsize*recordnum + SIZE_FLOAT64*2
    last = seg.firstword + SIZE_FLOAT64*seg.rsize*(recordnum+1)
    c = reinterpret(Float64, spk.daf.array[first:last], spk.daf.little)
    x = zeros(Float64, order)
    tc = 2.0 * frac/seg.intlen - 1.0
    x[1] = 1.0
    x[2] = tc
    twotc = tc + tc
    for i = 3:length(x)
        x[i] = twotc*x[i-1] - x[i-2]
    end
    reshape(c, (order, components)), x, seg.intlen, twotc
end

function position(c::Matrix, x::Vector)
    c'*x
end

function position(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc = getcoefficients(spk, seg, tdb, tdb2)
    position(c, x)
end

function velocity(c::Matrix, x::Vector, dt::Float64, twotc::Float64)
    n = length(x)
    t = zeros(n)
    t[2] = 1.0
    if n > 2
        t[3] = twotc + twotc
        for i = 4:n
            t[i] = twotc*t[i-1] - t[i-2] + x[i-1] + x[i-1]
        end
    end
    t *= 2.0
    t /= dt
    c'*t
end

function velocity(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc = getcoefficients(spk, seg, tdb, tdb2)
    velocity(c, x, dt, twotc)
end


function state(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc = getcoefficients(spk, seg, tdb, tdb2)
    r = position(c, x)
    v = velocity(c, x, dt, twotc)
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
