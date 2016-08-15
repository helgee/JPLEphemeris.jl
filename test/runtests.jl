using JPLEphemeris
using Base.Test
import Compat: ASCIIString

path = abspath(joinpath(dirname(@__FILE__), "..", "deps"))

include("basic.jl")

const AU = 0.149597870700000000e+09

const SPK_URL = Dict{Int, ASCIIString}(
    430 => "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp",
    405 => "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp",
)
const TEST_URL = Dict{Int, ASCIIString}(
    430 => "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/430/testpo.430",
    405 => "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/testpo.405",
)

function testephemeris(denum)
    ephem = SPK("$path/de$denum.bsp")
    if isfile("$path/testpo.$denum")
        lines = open(readlines, "$path/testpo.$denum")
    else
        error("Test file 'testpo.$denum' not found.")
    end

    start = findfirst(lines .== "EOT\n") + 1
    for l in lines[start:end]
        de, date, jd, target, center, index, value = split(l)
        target = parse(Int, target)
        center = parse(Int, center)
        index = parse(Int, index)
        value = float(value)
        jd = float(jd)

        if target in 14:15
            continue
        end

        try
            tr = teststate(ephem, jd, target)
            cr = teststate(ephem, jd, center)
            # From km/s to km/day
            tr[4:6] *= 86400
            cr[4:6] *= 86400
            # To AU and AU/day
            r = (tr - cr)/AU

            @test isapprox(r[index], value, atol=1e-13)
        # The test file for DE405 contains a wider range of dates than the SPK kernel provides.
        # Ignore those.
        catch err
            if isa(err, JPLEphemeris.OutOfRangeError)
                continue
            else
                rethrow(err)
            end
        end
    end
end

function teststate(kernel, tbd, target)
    if target == 3
        rv1 = state(kernel, 0, 3, tbd)
        rv2 = state(kernel, 3, 399, tbd)
        rv = rv1+rv2
    elseif target == 10
        rv1 = state(kernel, 0, 3, tbd)
        rv2 = state(kernel, 3, 301, tbd)
        rv = rv1+rv2
    elseif target == 11
        rv = state(kernel, 0, 10, tbd)
    elseif target == 12
        rv = zeros(6)
    elseif target == 13
        rv = state(kernel, 0, 3, tbd)
    else
        rv = state(kernel, 0, target, tbd)
    end
    return rv
end

# Run the JPL testsuite for every installed ephemeris.
@testset "Kernels" begin
    @testset for denum in ("430", "405")
        if !isdir(path)
            mkdir(path)
        end
        if !isfile("$path/de$denum.bsp")
            download(SPK_URL[denum], "$path/de$denum.bsp")
        end
        if !isfile("$path/testpo.$denum")
            download(TEST_URL[denum], "$path/testpo.$denum")
        end
        testephemeris(denum)
    end
end
