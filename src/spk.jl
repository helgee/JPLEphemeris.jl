import Base.position

const SECONDS_PER_DAY = 86400
const SIZE_FLOAT64 = sizeof(Float64)

type Segment
    name::ASCIIString
    firstsec::Float64
    lastsec::Float64
    firstdate::Float64
    lastdate::Float64
    target::Int32
    center::Int32
    frame::Int32
    repr::Int32
    firstaddr::Int32
    lastaddr::Int32
    firstword::Int32
    lastword::Int32
    intlen::Float64
    rsize::Int32
    n::Int32
end

function Segment(daf, name, record, little)
    firstsec, lastsec = reinterpret(Float64, record[1:16], little)
    target, center, frame, repr, firstaddr, lastaddr = reinterpret(Int32, record[17:end], little)
    if repr != 2
        error("Only Type 2 SPK files are supported.")
    end
    intlen, rsize, n = reinterpret(Float64, daf.array[lastaddr*SIZE_FLOAT64-3*SIZE_FLOAT64+1:lastaddr*SIZE_FLOAT64], little)
    Segment(
        name,
        firstsec,
        lastsec,
        firstsec/SECONDS_PER_DAY,
        lastsec/SECONDS_PER_DAY,
        target,
        center,
        frame,
        repr,
        firstaddr,
        lastaddr,
        firstaddr*SIZE_FLOAT64 - SIZE_FLOAT64 + 1,
        lastaddr*SIZE_FLOAT64 - SIZE_FLOAT64*4,
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
        seg = Segment(daf, name, summary, daf.little)
        if haskey(segments, seg.center)
            merge!(segments[seg.center], Dict(seg.target=>seg))
        else
            merge!(segments, Dict(seg.center=>Dict(seg.target=>seg)))
        end
    end
    SPK(daf, segments)
end

function checkdate(seg::Segment, tdb::Float64)
    if seg.firstdate > tdb > seg.lastdate
        error("The requested date '$tdb' is outside the range of this ephemeris.")
    end
end

function getcoefficients(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    checkdate(seg, tdb)
    components = 3
    order = (seg.rsize - 2) ÷ 3
    dt = seg.intlen/SECONDS_PER_DAY
    if tdb+tdb2 == seg.lastdate
        recordnum = seg.n - 1
        frac = dt
    else
        recordnum, frac = divrem(tdb - seg.firstdate, dt)
        recordnum = round(Int, recordnum)
    end
    # Drop the MID and RADIUS values
    first = seg.firstword + SIZE_FLOAT64*seg.rsize*recordnum + SIZE_FLOAT64*2
    last = seg.firstword + SIZE_FLOAT64*seg.rsize*(recordnum+1)
    c = reinterpret(Float64, spk.daf.array[first:last], spk.daf.little)
    x = zeros(Float64, order)
    tc = 2.0 * frac/dt - 1.0
    x[1] = 1.0
    x[2] = tc
    twotc = tc + tc
    for i = 3:length(x)
        x[i] = twotc*x[i-1] - x[i-2]
    end
    reshape(c, (order, components)), x, dt, twotc
end

function position(c::Matrix, x::Vector)
    c'*x
end

function position(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc = getcoefficients(spk, seg, tdb, tdb2)
    position(c, x)
end

function position(spk::SPK, center::Int, target::Int, tdb::Float64, tdb2::Float64=0.0)
    seg = spk.segments[center][target]
    position(spk, seg, tdb, tdb2)
end

function position(spk::SPK, target::Int, tdb::Float64, tdb2::Float64=0.0)
    seg = spk.segments[0][target]
    position(spk, seg, tdb, tdb2)
end

function velocity(c::Matrix, x::Vector, dt::Float64, twotc::Float64)
    t = zeros(x)
    t[2] = 1.0
    t[3] = twotc + twotc
    for i = 4:length(t)
        t[i] = twotc*t[i-1] - t[i-2] + x[i-1] + x[i-1]
    end
    t *= 2.0
    t /= dt
    c'*t
end

function velocity(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc = getcoefficients(spk, seg, tdb, tdb2)
    velocity(c, x, dt, twotc)
end

function velocity(spk::SPK, center::Int, target::Int, tdb::Float64, tdb2::Float64=0.0)
    seg = spk.segments[center][target]
    velocity(spk, seg, tdb, tdb2)
end

function velocity(spk::SPK, target::Int, tdb::Float64, tdb2::Float64=0.0)
    seg = spk.segments[0][target]
    velocity(spk, seg, tdb, tdb2)
end

function state(spk::SPK, seg::Segment, tdb::Float64, tdb2::Float64=0.0)
    c, x, dt, twotc = getcoefficients(spk, seg, tdb, tdb2)
    r = position(c, x)
    v = velocity(c, x, dt, twotc)
    [r;v]
end

function state(spk::SPK, center::Int, target::Int, tdb::Float64, tdb2::Float64=0.0)
    seg = spk.segments[center][target]
    state(spk, seg, tdb, tdb2)
end

function state(spk::SPK, target::Int, tdb::Float64, tdb2::Float64=0.0)
    seg = spk.segments[0][target]
    state(spk, seg, tdb, tdb2)
end
