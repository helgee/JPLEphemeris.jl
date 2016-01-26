const SECONDS_PER_DAY = 86400

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
end

type SPK
    daf::DAF
    segments::Dict{Int, Dict{Int, Segment}}
end

function Segment(name, record, little)
    firstsec, lastsec = reinterpret(Float64, record[1:16], little)
    target, center, frame, repr, firstaddr, lastaddr = reinterpret(Int32, record[17:end], little)
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
    )
end

function getdata(spk::SPK, seg::Segment)
    first = seg.firstaddr*8 - 7
    last = seg.lastaddr*8
    spk.daf.array[first:last]
end

function SPK(filename)
    daf = DAF(filename)
    segments = Dict{Int, Dict{Int, Segment}}()
    for (name, summary) in getsummaries(daf)
        seg = Segment(name, summary, daf.little)
        if haskey(segments, seg.center)
            merge!(segments[seg.center], Dict(seg.target=>seg))
        else
            merge!(segments, Dict(seg.center=>Dict(seg.target=>seg)))
        end
    end
    SPK(daf, segments)
end
