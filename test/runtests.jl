using JPLEphemeris
using Base.Test

const AU = 0.149597870700000000e+09

jd2000(jd) = jd - 2451544.5

function testephemeris(denum, verbose=false)
    println("Testing ephemeris DE$denum.")
    println(path)
    ephem = SPK("$path/de$denum.bsp")
    if isfile("$path/testpo.$denum")
        lines = open(readlines, "$path/testpo.$denum")
    else
        error("Test file 'testpo.$denum' not found.")
    end

    start = findfirst(lines .== "EOT\n") + 1
    for l in lines[start:end]
        de, date, jd, target, center, index, value = split(l)
        jd = jd2000(float(jd))
        target = parse(Int, target)
        center = parse(Int, center)
        index = parse(Int, index)
        value = float(value)

        if target in 14:15
            continue
        end

        tr = teststate(ephem, jd, target)
        cr = teststate(ephem, jd, center)
        r = (tr - cr)/AU

        if verbose
            println("Date: $date")
            println("JD2000: $jd")
            println("Target: $target")
            println("Center: $center")
            println("Orginal value: $value")
            println("Computed value: $(r[index])")
            println("===========================================")
        end

        #= @test_approx_eq r[index] value =#
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

path = abspath(joinpath(splitdir(@__FILE__)[1], "..", "deps"))

verbose = false
if ~isempty(ARGS) && ((ARGS[1] == "-v") || (ARGS[1] == "--verbose"))
    verbose = true
end

# Run the JPL testsuite for every installed ephemeris.
for denum in (430, 405)
    testephemeris(denum, verbose)
end

