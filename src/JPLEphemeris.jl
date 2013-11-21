module JPLEphemeris

using HDF5, JLD

export PLANETS

type Ephemeris
    startepoch::Float64
    finalepoch::Float64
    file::String
    groups::Array{Float64,2}
end

const PLANETS = [
    "Mercury",
    "Venus",
    "Earth-Moon Barycenter",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "Pluto",
    "Moon (geocentric)",
    "Sun",
    "Nutations",
    "Librations",
    ]

#function readascii(header, datafiles, hdf5file)
function readascii(header, outfile)
    startepoch = 0.0
    finalepoch = 0.0
    constantnames = String[]
    constantvalues = Float64[]
    l = open(readlines, header)
    for i = 1:length(l)
        if beginswith(l[i], "GROUP   1030")
            startepoch, finalepoch = split(l[i+2])[1:2]
        elseif beginswith(l[i], "GROUP   1040")
            n = int(l[i+2])
            firstline = i+3
            lastline = i+3+div(n,10)
            for j = firstline:lastline
                append!(constantnames, split(l[j]))
            end
        elseif beginswith(l[i], "GROUP   1041")
            n = int(l[i+2])
            firstline = i+3
            lastline = i+3+div(n,3)
            for j = firstline:lastline
                append!(constantvalues, float(split(replace(l[j],"D","e"))))
            end
        end
    end
    groupreg = r"^\s+[0-9]\s+[0-9]+"
    constants = Dict(constantnames, constantvalues)
    file = jldopen(outfile, "w")
    @write file constants
    close(file)
end

function poly(a::Vector, x::Float64)
    y = 0.0
    for i = 1:length(a)
        y += a[i]*x^i
    end
    return y
end

end
