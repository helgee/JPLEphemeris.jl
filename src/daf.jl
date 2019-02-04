using Mmap

struct DAF
    filename::String
    array::Vector{UInt8}
    little::Bool
    id::String
    nd::Int32
    ni::Int32
    name::String
    first::Int32
    last::Int32
    ss::Int32
    nc::Int32
end

const RECORD_LENGTH = 1024

function getrecord(array, n)
    array[1+RECORD_LENGTH*(n-1):RECORD_LENGTH*n]
end

getrecord(daf::DAF, n) = getrecord(daf.array, n)

function DAF(filename)
    if !isfile(filename)
        throw(ArgumentError("'$filename' does not exist."))
    end
    array = Mmap.mmap(filename, Vector{UInt8})
    fr = readfilerecord(getrecord(array, 1))
    DAF(filename, array, fr...)
end

@inline function reinterpret_getindex(::Type{T}, b::Vector{UInt8}, i, le::Bool) where T
    # Read in host endian
    @boundscheck if i + sizeof(T) - 1 > length(b)
        throw(BoundsError(b, i + sizeof(T) - 1))
    end
    v = unsafe_load(Ptr{T}(pointer(b, i)))
    le ? htol(v) : hton(v)
end

@inline function reinterpret_getindex(::Type{T}, b::Vector{UInt8},
                                      idxs::NTuple{N,Int}, le::Bool) where {N, T}
    @inbounds ntuple(i->reinterpret_getindex(T, b, idxs[i], le), Val(N))
end

function readint(record, address, littleendian=true)
    reinterpret_getindex(Int32, record, address + 1, littleendian)
end

function readascii(record, address, len)
    rstrip(String(record[address+1:address+len]))
end

function islittleendian(record, legacy)
    if legacy
        nd = readint(record, 8)
        if nd == 2
            return true
        end
        nd = readint(record, 8, false)
        if nd == 2
            return false
        else
            error("Endianess could not be detected.")
        end
    end

    endianess = readascii(record, 88, 8)
    endianess == "LTL-IEEE"
end

function readfilerecord(record)
    ftpstr = b"FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP"
    id = readascii(record,0,8)
    legacy = id == "NAIF/DAF"
    little = islittleendian(record, legacy)
    nd = readint(record, 8, little)
    ni = readint(record, 12, little)
    name = readascii(record, 16, 60)
    first = readint(record, 76, little)
    last = readint(record, 80, little)
    ss = nd + div(ni+1, 2)
    nc = 8*ss
    if ~legacy
        if ftpstr != record[699+1:699+28]
            error("This DAF file is damaged.")
        end
    end
    little, id, nd, ni, name, first, last, ss, nc
end

function summaryheader(record, little)
    next, nsum = reinterpret_getindex(Float64, record, (1, 17), little)
    round(Int32, next), round(Int32, nsum)
end

function addsummaries!(summaries, record, names, nsum, ss, nc)
    for i = 1:nsum
        push!(summaries, (
            rstrip(String(names[1+(i-1)*nc:i*nc])),
            record[25+(i-1)*ss*8:25+i*ss*8]
        ))
    end
end

function getsummaries(daf::DAF)
    summaries = Tuple{String,Vector{UInt8}}[]
    record = getrecord(daf, daf.first)
    names = getrecord(daf, daf.first+1)
    next, nsum = summaryheader(record, daf.little)
    addsummaries!(summaries, record, names, nsum, daf.ss, daf.nc)
    while next != 0
        record = getrecord(daf, next)
        names = getrecord(daf, next+1)
        next, nsum = summaryheader(record, daf.little)
        addsummaries!(summaries, record, names, nsum, daf.ss, daf.nc)
    end
    summaries
end
