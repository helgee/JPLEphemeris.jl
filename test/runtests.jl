if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
using JPLEphemeris
import Compat: ASCIIString

const AU = 0.149597870700000000e+09

const SPK_URL = Dict{Int, ASCIIString}(
    430 => "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp",
    405 => "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp",
)
const TEST_URL = Dict{Int, ASCIIString}(
    430 => "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/430/testpo.430",
    405 => "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/testpo.405",
)

path = abspath(joinpath(dirname(@__FILE__), "..", "deps"))

for denum in (430, 405)
    if !isdir(path)
        mkdir(path)
    end
    if !isfile("$path/de$denum.bsp")
        download(SPK_URL[denum], "$path/de$denum.bsp")
    end
    if !isfile("$path/testpo.$denum")
        download(TEST_URL[denum], "$path/testpo.$denum")
    end
end

include("basic.jl")
include("kernels.jl")
