# Reference value from CSPICE
rvm = [4.250906022073639e7,2.3501057648129586e7,8.158467467032234e6,-34.160008371029825,37.844059275357594,23.756128199757867]

de430segments = [
    "SOLAR SYSTEM BARYCENTER (0) => MERCURY BARYCENTER (1)",
    "SOLAR SYSTEM BARYCENTER (0) => VENUS BARYCENTER (2)",
    "SOLAR SYSTEM BARYCENTER (0) => EARTH-MOON BARYCENTER (3)",
    "SOLAR SYSTEM BARYCENTER (0) => MARS BARYCENTER (4)",
    "SOLAR SYSTEM BARYCENTER (0) => JUPITER BARYCENTER (5)",
    "SOLAR SYSTEM BARYCENTER (0) => SATURN BARYCENTER (6)",
    "SOLAR SYSTEM BARYCENTER (0) => URANUS BARYCENTER (7)",
    "SOLAR SYSTEM BARYCENTER (0) => NEPTUNE BARYCENTER (8)",
    "SOLAR SYSTEM BARYCENTER (0) => PLUTO BARYCENTER (9)",
    "SOLAR SYSTEM BARYCENTER (0) => SUN (10)",
    "MERCURY BARYCENTER (1) => MERCURY (199)",
    "VENUS BARYCENTER (2) => VENUS (299)",
    "EARTH-MOON BARYCENTER (3) => MOON (301)",
    "EARTH-MOON BARYCENTER (3) => EARTH (399)",
]
jd = Dates.datetime2julian(DateTime(2016,1,1,0,0,0))
spk = SPK("$path/de430.bsp")

@testset "API" begin
    @testset for (a, b) in zip(JPLEphemeris.list_segments(spk), de430segments)
        @test a == b
    end
    @testset for (func, range) in zip((position, velocity, state), (1:3, 4:6, 1:6))
        res = func(spk, "mercury barycenter", jd)
        @test res ≈ rvm[range]
        res = func(spk, 1, jd)
        @test res ≈ rvm[range]
        res = func(spk, "ssb", "mercury barycenter", jd)
        @test res ≈ rvm[range]
        res = func(spk, 0, "mercury barycenter", jd)
        @test res ≈ rvm[range]
        res = func(spk, "ssb", 1, jd)
        @test res ≈ rvm[range]
        res = func(spk, 0, 1, jd)
        @test res ≈ rvm[range]
        res = func.(spk, 0, 1, [jd; jd])
        @test all(map(x -> x ≈ rvm[range], res))
        res = func(spk, "mercury barycenter", jd, 0.0)
        @test res ≈ rvm[range]
        res = func(spk, 1, jd, 0.0)
        @test res ≈ rvm[range]
        res = func(spk, "ssb", "mercury barycenter", jd, 0.0)
        @test res ≈ rvm[range]
        res = func(spk, 0, "mercury barycenter", jd, 0.0)
        @test res ≈ rvm[range]
        res = func(spk, "ssb", 1, jd, 0.0)
        @test res ≈ rvm[range]
        res = func(spk, 0, 1, jd, 0.0)
        @test res ≈ rvm[range]
        res = func.(spk, 0, 1, [jd; jd], [0.0, 0.0])
        @test all(map(x -> x ≈ rvm[range], res))
    end
end
