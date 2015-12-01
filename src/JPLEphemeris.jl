module JPLEphemeris

using HDF5, JLD

import Base.position

export Ephemeris, STANDARD_EPHEMERIS
export position, velocity, state, close, getephem, rmephem

import JLD.close

const STANDARD_EPHEMERIS = "430"
const PATH = "$(Pkg.dir())/JPLEphemeris/deps"


type Ephemeris
    id::Int
    fid::JLD.JldFile
    startepoch::Float64
    finalepoch::Float64
    dtable::Vector{Vector{Float64}}
    intervall::Vector{Float64}
    cache::Dict{AbstractString, Matrix{Float64}}
    constants::Dict{AbstractString, Float64}

    function Ephemeris(file::AbstractString)
        fid = jldopen(file, "r")
        startepoch = read(fid, "startepoch")
        finalepoch = read(fid, "finalepoch")
        dtable = read(fid, "dtable")
        intervall = read(fid, "intervall")
        cache = Dict{AbstractString, Matrix{Float64}}()
        sizehint!(cache, mapreduce(length, +, dtable))
        constants = read(fid, "constants")
        id = constants["DENUM"]
        merge!(constants, Dict("earthshare"=>1.0/(1.0+constants["EMRAT"]),
        "moonshare"=>constants["EMRAT"]/(1.0+constants["EMRAT"])))
        return new(id, fid, startepoch, finalepoch, dtable, intervall, cache, constants)
    end
end

function Ephemeris()
    return Ephemeris("$PATH/de$STANDARD_EPHEMERIS.jld")
end

function Ephemeris(denum::Integer)
    if isfile("$PATH/de$denum.jld")
        return Ephemeris("$PATH/de$denum.jld")
    else
        error("Ephemeris file '$PATH/de$denum.jld' not found.")
    end
end

const planets = Dict(
    "mercury"=>1,
    "venus"=>2,
    "earthmoon"=>3,
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
)

function position(ephem::Ephemeris, body::AbstractString, date::Float64)
    i = planets[body]
    checkdate(ephem, date)
    c = coefficients(ephem, i, date)
    return position(c...)
end

function position(ephem::Ephemeris, body::AbstractString, date::AbstractVector{Float64})
    n = body == "nutations" ? 2 : 3
    p = zeros(n, length(date))
    for (i,d) in enumerate(date)
        p[:,i] = position(ephem, body, d)
    end
    return p
end

function velocity(ephem::Ephemeris, body::AbstractString, date::Float64)
    i = planets[body]
    checkdate(ephem, date)
    c = coefficients(ephem, i, date)
    return velocity(c...)
end

function velocity(ephem::Ephemeris, body::AbstractString, date::AbstractVector{Float64})
    n = body == "nutations" ? 2 : 3
    v = zeros(n, length(date))
    for (i,d) in enumerate(date)
        v[:,i] = velocity(ephem, body, d)
    end
    return v
end

function state(ephem::Ephemeris, body::AbstractString, date::Float64)
    nbody = planets[body]
    checkdate(ephem, date)
    c = coefficients(ephem, nbody, date)
    return [position(c...); velocity(c...)]
end

function state(ephem::Ephemeris, body::AbstractString, date::AbstractVector{Float64})
    n = body == "nutations" ? 4 : 6
    s = zeros(n, length(date))
    for (i,d) in enumerate(date)
        s[:,i] = state(ephem, body, d)
    end
    return s
end

function checkdate(ephem::Ephemeris, date::Float64)
    if (date < ephem.startepoch) | (date > ephem.finalepoch)
        error("The date must be between $(ephem.startepoch) and
        $(ephem.finalepoch).")
    end
end

function coefficients(ephem::Ephemeris, nbody::Int64, date::Float64)
    dt = ephem.intervall[nbody]
    if date == ephem.finalepoch
        frac = dt
        d = ephem.dtable[nbody][end]
    else
        index, frac = divrem(date-ephem.startepoch, dt)
        index = round(Int,(index) + 1)
        d = ephem.dtable[nbody][index] 
    end

    # Check if the requested coeffcients are already cached,
    # if not retrieve them from the JLD file.
    key = "$nbody-$d"
    if ~haskey(ephem.cache, key)
        merge!(ephem.cache, Dict(key=>read(ephem.fid, key)))
    end
    c = ephem.cache[key]
    x = zeros(size(c)[1])

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

function position(c::Matrix{Float64}, x::Vector{Float64},
    dt::Float64, twotc::Float64)
    return c'*x
end

function velocity(c::Matrix{Float64}, x::Vector{Float64},
    dt::Float64, twotc::Float64)
    t = zeros(length(x))
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

include("util.jl")

end
