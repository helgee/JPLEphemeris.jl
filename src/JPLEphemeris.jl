module JPLEphemeris

using HDF5, JLD

export Ephemeris
export position, velocity, state

type Ephemeris
    fid::JLD.JldFile
    startepoch::Float64
    finalepoch::Float64
    dtable::Vector{Vector{Float64}}
    intervall::Vector{Float64}
    cache::Dict{String, Matrix{Float64}}
    constants::Dict{String, Float64}

    function Ephemeris(file::String)
        fid = jldopen(file, "r")
        startepoch = read(fid, "startepoch")
        finalepoch = read(fid, "finalepoch")
        dtable = read(fid, "dtable")
        intervall = read(fid, "intervall")
        cache = Dict{String, Matrix{Float64}}()
        sizehint(cache, maximum([length(dtable[i]) for i = 1:length(dtable)]))
        constants = read(fid, "constants")
        return new(fid, startepoch, finalepoch, dtable, intervall, cache, constants)
    end
end

const planets = [
    "mercury"=>1,
    "venus"=>2,
    "earth"=>3,
    "mars"=>4,
    "jupiter"=>5,
    "saturn"=>6,
    "uranus"=>7,
    "neptune"=>8,
    "pluto"=>9,
    "moon"=>10,
    "sun"=>11,
    "nutations"=>12,
    "librations"=>13,
]

function position(ephem::Ephemeris, body::String, date::Float64)
    i = planets[body]
    if (date < ephem.startepoch) | (date > ephem.finalepoch)
        error("The date must be between $(ephem.startepoch) and
        $(ephem.finalepoch).")
    end
    c = coefficients(ephem, i, date)
    return pos(c...)
end

function position(ephem::Ephemeris, body::String, date::Vector{Float64})
    n = body == "nutations" ? 2 : 3
    p = zeros(Float64, n, length(date))
    for (i,d) in enumerate(date)
        p[:,i] = position(ephem, body, d)
    end
    return p
end

function velocity(ephem::Ephemeris, body::String, date::Float64)
    i = planets[body]
    if (date < ephem.startepoch) | (date > ephem.finalepoch)
        error("The date must be between $(ephem.startepoch) and
        $(ephem.finalepoch).")
    end
    c = coefficients(ephem, i, date)
    return vel(c...)
end

function velocity(ephem::Ephemeris, body::String, date::Vector{Float64})
    n = body == "nutations" ? 2 : 3
    v = zeros(Float64, n, length(date))
    for (i,d) in enumerate(date)
        v[:,i] = velocity(ephem, body, d)
    end
    return v
end

function state(ephem::Ephemeris, body::String, date::Float64)
    nbody = planets[body]
    if (date < ephem.startepoch) | (date > ephem.finalepoch)
        error("The date must be between $(ephem.startepoch) and
        $(ephem.finalepoch).")
    end
    c = coefficients(ephem, nbody, date)
    return [pos(c...), vel(c...)]
end

function state(ephem::Ephemeris, body::String, date::Vector{Float64})
    n = body == "nutations" ? 4 : 6
    s = zeros(Float64, n, length(date))
    for (i,d) in enumerate(date)
        s[:,i] = state(ephem, body, d)
    end
    return s
end

function coefficients(ephem::Ephemeris, nbody::Int64, date::Float64)
    dt = ephem.intervall[nbody]
    if date == ephem.finalepoch
        frac = dt
        d = ephem.dtable[nbody][end]
    else
        index, frac = divrem(date-ephem.startepoch, dt)
        index = int(index) + 1
        d = ephem.dtable[nbody][index] 
    end

    # Check if the requested coeffcients are already cached,
    # if not retrieve them from the HDF5 file.
    key = "$nbody-$d"
    if ~haskey(ephem.cache, key)
        merge!(ephem.cache, [key=>read(ephem.fid, key)])
    end
    c = ephem.cache[key]
    x = similar(c, size(c)[1])

    # Normalized Chebyshev time
    tc = 2.0 * frac/dt - 1.0
    x[1] = 1
    x[2] = tc
    twotc = tc + tc
    for i = 3:length(x)
        x[i] = twotc*x[i-1] - x[i-2]
    end
    return c, x, dt, twotc
end

function pos(c, x, dt, twotc)
    return c'*x
end

function vel(c, x, dt, twotc)
    t = similar(x, length(x))
    t[1] = 0.0
    t[2] = 1.0
    t[3] = twotc + twotc
    for i = 4:length(t)
        t[i] = twotc*t[i-1] - t[i-2] + x[i-1] + x[i-1]
    end
    t *= 2.0
    t /= dt
    return c'*t
end

function close(ephem::Ephemeris)
    close(ephem.fid)
end

include("parser.jl")

end
