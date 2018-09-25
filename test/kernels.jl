function testephemeris(denum)
    ephem = SPK("$path/de$denum.bsp")
    if isfile("$path/testpo.$denum")
        lines = open(readlines, "$path/testpo.$denum")
    else
        error("Test file 'testpo.$denum' not found.")
    end

    start = findfirst(map(x -> startswith(x, "EOT"), lines)) + 1
    for l in lines[start:end]
        de, date, jd, target, center, index, value = split(l)
        target = parse(Int, target)
        center = parse(Int, center)
        index = parse(Int, index)
        value = parse(Float64, value)
        jd = parse(Float64, jd)

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

            passed = isapprox(r[index], value, atol=1e-8)
            if !passed
                @show jd
                @show tr
                @show cr
                @show r
                @show target
                @show center
                @show index
                @show r[index]
                @show value
                @show abs(r[index] - value)
            end
            @test passed
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
        rv = rv1 .+ rv2
    elseif target == 10
        rv1 = state(kernel, 0, 3, tbd)
        rv2 = state(kernel, 3, 301, tbd)
        rv = rv1 .+ rv2
    elseif target == 11
        rv = state(kernel, 0, 10, tbd)
    elseif target == 12
        rv = zeros(6)
    elseif target == 13
        rv = state(kernel, 0, 3, tbd)
    else
        rv = state(kernel, 0, target, tbd)
    end
    return vcat(rv...)
end

# Run the JPL testsuite for every installed ephemeris.
@testset "Kernels" begin
    @testset for denum in (430, 405)
        testephemeris(denum)
    end
end
