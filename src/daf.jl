type DAF
    filename::AbstractString
    array::Vector{UInt8}
    LOCIDW::ASCIIString
    ND::Int32
    NI::Int32
    LOCIFN::ASCIIString
    FWARD::Int32
    BWARD::Int32
    FREE::Int32
    LOCFMT::ASCIIString
    SS::Int32
    NC::Int32
end

const RECORD_LENGTH = 1024

function getrecord(array, n)
    array[1+RECORD_LENGTH*(n-1):RECORD_LENGTH*n]
end

getrecord(daf::DAF, n) = getrecord(daf.array, n)

function DAF(filename)
    array = Mmap.mmap(filename, Vector{UInt8})
    fr = readfilerecord(getrecord(array, 1))
    DAF(
        filename,
        array,
        fr["LOCIDW"],
        fr["ND"],
        fr["NI"],
        fr["LOCIFN"], 
        fr["FWARD"],
        fr["BWARD"],
        fr["FREE"],
        fr["LOCFMT"],
        fr["SS"],
        fr["NC"],
    )
end

address(start, len) = start+1:start+len

function readfilerecord(record)
    FTPSTR = b"FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP"
    LOCIDW = strip(ascii(record[address(0,8)]))
    ND = reinterpret(Int32, record[address(8,4)])[1]
    NI = reinterpret(Int32, record[address(12,4)])[1]
    LOCIFN = strip(ascii(record[address(16,60)]))
    FWARD = reinterpret(Int32, record[address(76,4)])[1]
    BWARD = reinterpret(Int32, record[address(80,4)])[1]
    FREE = reinterpret(Int32, record[address(84,4)])[1]
    LOCFMT = strip(ascii(record[address(88,8)]))
    SS = ND + div(NI+1,2)
    NC = 8*SS
    if FTPSTR != record[address(699,28)]
        error("This DAF file is damaged.")
    end
    return Dict(
        "LOCIDW" => LOCIDW,
        "ND" => ND,
        "NI" => NI,
        "LOCIFN" => LOCIFN,
        "FWARD" => FWARD,
        "BWARD" => BWARD,
        "FREE" => FREE,
        "LOCFMT" => LOCFMT,
        "SS" => SS,
        "NC" => NC,
    )
end

function getsummaries(daf::DAF)
    summaries, nextsummary = getsummaries(getrecord(daf, daf.FWARD), getrecord(daf, daf.FWARD+1), daf.SS, daf.NC)
    while nextsummary != 0
        s, nextsummary = getsummaries(getrecord(daf, nextsummary), getrecord(daf, nextsummary+1), daf.SS, daf.NC)
        append!(summaries, s)
    end
    return summaries
end

function getsummaries(summaryrecord, namesrecord, summarylength, namelength)
    nextsummary, _, nsum = round(Int32, reinterpret(Float64, summaryrecord[1:24]))
    summaries = Tuple{ASCIIString,Vector{UInt8}}[]
    for i = 1:nsum
        push!(summaries,(
            rstrip(ascii(namesrecord[1+(i-1)*namelength:i*namelength])),
            summaryrecord[25+(i-1)*summarylength*8:25+i*summarylength*8]
        ))
    end
    return summaries, nextsummary
end
