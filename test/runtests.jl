using JPLEphemeris
using Base.Test

import JPLEphemeris.state

function testephemeris(ephem::Ephemeris)
    denum = int(ephem.constants["DENUM"])
    if isfile("$path/testpo.$denum")
        lines = open(readlines, "$path/testpo.$denum")
    else
        error("Test file 'testpo.$denum' not found.")
    end

    for l in lines
        if ismatch(r"^[0-9]", l)
            de, date, jd, target, center, index, value = [f(v) for (f,v) in
            zip((string, string, float, int, int, int, float), split(l))]

            if target in 14:15
                r = state(ephem, jd, target)
            else
                tr = state(ephem, jd, target)
                cr = state(ephem, jd, center)
                r = (tr - cr)/ephem.constants["AU"]
            end

            #println("Date: $date")
            #println("Julian day: $jd")
            #println("Target: $target")
            #println("Center: $center")
            #println("Orginal value: $value")
            #println("Computed value: $(r[index])")
            #println("===========================================")

            if target == 15 && index == 3
                delta = (r[index] - value)/(0.23*(jd - 2451545.0))
                @test delta <= 1e-13
            elseif target == 15 && index == 6
                delta = (r[index] - value)*0.01/(1.0 + (jd - 2451545.0)/365.25)
                @test delta <= 1e-13
            else
                @test isapprox(r[index], value)
            end
        end
    end
end

function state(ephem::Ephemeris, date::Float64, target::Int64)
    s(body::String) = state(ephem, body, date)

    planets = [1=>"mercury", 2=>"venus", 4=>"mars", 5=>"jupiter",
    6=>"saturn", 7=>"uranus", 8=>"neptune", 9=>"pluto", 11=>"sun",
    14=>"nutations", 15=>"librations"]

    if target == 3
        return s("earthmoon") - s("moon") * ephem.constants["earthshare"]
    elseif target == 10
        return s("earthmoon") + s("moon") * ephem.constants["moonshare"]
    elseif target == 12
        return zeros(6)
    elseif target == 13
        return s("earthmoon")
    else
        return s(planets[target])
    end
end

path = "$(Pkg.dir())/JPLEphemeris/deps"
files = readdir(path)
length(files) < 2 && error("No ephemeris files installed.")

# Run the JPL testsuite for every installed ephemeris.
for f in files
    file, ext = splitext(f)
    if ext == ".jld"
        eph = Ephemeris("$path/$f")
        testephemeris(eph)
    end
end

