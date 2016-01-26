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

function Segment(name, record)
    firstsec, lastsec = reinterpret(Float64, record[1:16])
    target, center, frame, repr, firstaddr, lastaddr = reinterpret(Int32, record[17:end])
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

type SPK
    daf::DAF
    segments::Vector{Segment}
end

function SPK(filename)
    daf = DAF(filename)
    segments = [Segment(s...) for s in getsummaries(daf)]
    SPK(daf, segments)
end
