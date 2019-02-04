import AstroBase: TDBEpoch, mercury_barycenter, ssb, earth, moon,
    earth_barycenter, mercury

# Reference value from CSPICE
r_ref = [4.250906022073639e7, 2.3501057648129586e7, 8.158467467032234e6]
v_ref = [-34.160008371029825, 37.844059275357594, 23.756128199757867]
rv_ref = (r_ref, v_ref)

de430segments = [
    "Solar System Barycenter (0) => Mercury Barycenter (1)",
    "Solar System Barycenter (0) => Venus Barycenter (2)",
    "Solar System Barycenter (0) => Earth Barycenter (3)",
    "Solar System Barycenter (0) => Mars Barycenter (4)",
    "Solar System Barycenter (0) => Jupiter Barycenter (5)",
    "Solar System Barycenter (0) => Saturn Barycenter (6)",
    "Solar System Barycenter (0) => Uranus Barycenter (7)",
    "Solar System Barycenter (0) => Neptune Barycenter (8)",
    "Solar System Barycenter (0) => Pluto Barycenter (9)",
    "Solar System Barycenter (0) => Sun (10)",
    "Mercury Barycenter (1) => Mercury (199)",
    "Venus Barycenter (2) => Venus (299)",
    "Earth Barycenter (3) => Luna (301)",
    "Earth Barycenter (3) => Earth (399)",
]
ep = TDBEpoch(2016, 1, 1)
spk = SPK("$path/de430.bsp")

@testset "API" begin
    @testset "Segments" begin
        @test JPLEphemeris.findsegment(spk.segments, 0, 3) == (spk.segments[0][3], 1.0)
        @test JPLEphemeris.findsegment(spk.segments, 3, 0) == (spk.segments[0][3], -1.0)
    end
    @testset for (a, b) in zip(JPLEphemeris.list_segments(spk), de430segments)
        @test a == b
    end
    @testset for (func, ref) in zip((position, velocity, position_velocity), (r_ref, v_ref, rv_ref))
        res = func(spk, ep, ssb, mercury_barycenter)
        @test all(res .≈ ref)
        res = func(spk, ep, mercury_barycenter, ssb)
        @test all(res .≈ -1 .* ref)

        res = func(spk, ep, earth, mercury)
        exp = func(spk, ep, earth, earth_barycenter) .+
            func(spk, ep, earth_barycenter, ssb) .+
            func(spk, ep, ssb, mercury_barycenter) .+
            func(spk, ep, mercury_barycenter, mercury)
        @test all(res .≈ exp)
    end
end
