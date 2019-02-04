using AstroBase: TDBEpoch, from_naifid, ssb, sun, luna, earth_barycenter, earth

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
        jd = TDBEpoch(parse(Float64, jd), origin=:julian)

        if target in 14:15
            continue
        end

        try
            r = teststate(ephem, jd, center, target)

            passed = isapprox(r[index], value, atol=1e-8)
            if !passed
                @show jd
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

function id_to_body(id)
    id == 3 && return earth
    id == 10 && return luna
    id == 11 && return sun
    id == 12 && return ssb
    id == 13 && return earth_barycenter

    from_naifid(id)
end

function teststate(kernel, tdb, origin, target)
    origin = id_to_body(origin)
    target = id_to_body(target)
    rv = vcat(position_velocity(kernel, tdb, origin, target)...)
    # From km/s to km/day
    rv[4:6] .*= 86400.0
    # From km to AU
    rv./AU
end

# Run the JPL testsuite for every installed ephemeris.
@testset "Kernels" begin
    @testset for denum in (430, 405)
        testephemeris(denum)
    end
end
