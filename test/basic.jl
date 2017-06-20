import AstroDynBase: TDBEpoch, MercuryBarycenter, SSB, Earth, Moon,
    EarthBarycenter, Mercury

# Reference value from CSPICE
r_ref = [4.250906022073639e7, 2.3501057648129586e7, 8.158467467032234e6]
v_ref = [-34.160008371029825, 37.844059275357594, 23.756128199757867]
rv_ref = [r_ref; v_ref]

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
ep = TDBEpoch(jd)
spk = SPK("$path/de430.bsp")

@testset "API" begin
    @testset "Path" begin
        @test JPLEphemeris.findpath(EarthBarycenter, SSB) == [EarthBarycenter, SSB]
        @test JPLEphemeris.findpath(SSB, EarthBarycenter) == [SSB, EarthBarycenter]
        @test JPLEphemeris.findpath(Earth, SSB) == [Earth, EarthBarycenter, SSB]
        @test JPLEphemeris.findpath(SSB, Earth) == [SSB, EarthBarycenter, Earth]
        @test JPLEphemeris.findpath(Earth, Mercury) == [Earth, EarthBarycenter, SSB, MercuryBarycenter, Mercury]
        @test JPLEphemeris.findpath(Earth, Moon) == [Earth, EarthBarycenter, Moon]
        @test JPLEphemeris.findpath(EarthBarycenter, MercuryBarycenter) == [EarthBarycenter, SSB, MercuryBarycenter]
        @test JPLEphemeris.findpath(Earth, MercuryBarycenter) == [Earth, EarthBarycenter, SSB, MercuryBarycenter]
        @test JPLEphemeris.findpath(EarthBarycenter, Mercury) == [EarthBarycenter, SSB, MercuryBarycenter, Mercury]
    end
    @testset "Segments" begin
        @test JPLEphemeris.findsegment(spk.segments, 0, 3) == (spk.segments[0][3], 1.0)
        @test JPLEphemeris.findsegment(spk.segments, 3, 0) == (spk.segments[0][3], -1.0)
    end
    @testset for (a, b) in zip(JPLEphemeris.list_segments(spk), de430segments)
        @test a == b
    end
    @testset for (func, ref) in zip((position, velocity, state), (r_ref, v_ref, rv_ref))
        res = func(spk, "mercury barycenter", jd)
        @test all(res .≈ ref)
        res = func(spk, 1, jd)
        @test all(res .≈ ref)
        res = func(spk, "ssb", "mercury barycenter", jd)
        @test all(res .≈ ref)
        res = func(spk, 0, "mercury barycenter", jd)
        @test all(res .≈ ref)
        res = func(spk, "ssb", 1, jd)
        @test all(res .≈ ref)
        res = func(spk, 0, 1, jd)
        @test all(res .≈ ref)
        res = func.(spk, 0, 1, [jd; jd])
        @test all(all.(map(x->x .≈ ref, res)))
        res = func(spk, "mercury barycenter", jd, 0.0)
        @test all(res .≈ ref)
        res = func(spk, 1, jd, 0.0)
        @test all(res .≈ ref)
        res = func(spk, "ssb", "mercury barycenter", jd, 0.0)
        @test all(res .≈ ref)
        res = func(spk, 0, "mercury barycenter", jd, 0.0)
        @test all(res .≈ ref)
        res = func(spk, "ssb", 1, jd, 0.0)
        @test all(res .≈ ref)
        res = func(spk, 0, 1, jd, 0.0)
        @test all(res .≈ ref)
        res = func.(spk, 0, 1, [jd; jd], [0.0, 0.0])
        @test all(all.(map(x->x .≈ ref, res)))

        res = func(spk, ep, SSB, MercuryBarycenter)
        @test all(res .≈ ref)
        res = func(spk, ep, MercuryBarycenter, SSB)
        @test all(res .≈ -1 .* ref)

        res = func(spk, ep, Earth, Mercury)
        exp = func(spk, ep, Earth, EarthBarycenter) .+
            func(spk, ep, EarthBarycenter, SSB) .+
            func(spk, ep, SSB, MercuryBarycenter) .+
            func(spk, ep, MercuryBarycenter, Mercury)
        @test all(res .≈ exp)
    end
end
