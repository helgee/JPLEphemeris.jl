using Base.Test
using JPLEphemeris

# Reference value from CSPICE
rvm = [4.250906022073639e7,2.3501057648129586e7,8.158467467032234e6,-34.160008371029825,37.844059275357594,23.756128199757867]
mrvm = [rvm';rvm']

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

spk = SPK("$path/de430.bsp")
@test JPLEphemeris.list_segments(spk) == de430segments
jd = Dates.datetime2julian(DateTime(2016,1,1,0,0,0))
pos = position(spk, "mercury barycenter", jd)
@test pos ≈ rvm[1:3]
pos = position(spk, 1, jd)
@test pos ≈ rvm[1:3]
pos = position(spk, "ssb", "mercury barycenter", jd)
@test pos ≈ rvm[1:3]
pos = position(spk, 0, "mercury barycenter", jd)
@test pos ≈ rvm[1:3]
pos = position(spk, "ssb", 1, jd)
@test pos ≈ rvm[1:3]
pos = position(spk, 0, 1, jd)
@test pos ≈ rvm[1:3]
pos = position(spk, 0, 1, [jd; jd])
@test pos ≈ mrvm[:,1:3]
pos = position(spk, "mercury barycenter", jd, 0.0)
@test pos ≈ rvm[1:3]
pos = position(spk, 1, jd, 0.0)
@test pos ≈ rvm[1:3]
pos = position(spk, "ssb", "mercury barycenter", jd, 0.0)
@test pos ≈ rvm[1:3]
pos = position(spk, 0, "mercury barycenter", jd, 0.0)
@test pos ≈ rvm[1:3]
pos = position(spk, "ssb", 1, jd, 0.0)
@test pos ≈ rvm[1:3]
pos = position(spk, 0, 1, jd, 0.0)
@test pos ≈ rvm[1:3]
pos = position(spk, 0, 1, [jd; jd], [0.0, 0.0])
@test pos ≈ mrvm[:,1:3]

vel = velocity(spk, "mercury barycenter", jd)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 1, jd)
@test vel ≈ rvm[4:6]
vel = velocity(spk, "ssb", "mercury barycenter", jd)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 0, "mercury barycenter", jd)
@test vel ≈ rvm[4:6]
vel = velocity(spk, "ssb", 1, jd)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 0, 1, jd)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 0, 1, [jd; jd])
@test vel ≈ mrvm[:,4:6]
vel = velocity(spk, "mercury barycenter", jd, 0.0)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 1, jd, 0.0)
@test vel ≈ rvm[4:6]
vel = velocity(spk, "ssb", "mercury barycenter", jd, 0.0)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 0, "mercury barycenter", jd, 0.0)
@test vel ≈ rvm[4:6]
vel = velocity(spk, "ssb", 1, jd, 0.0)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 0, 1, jd, 0.0)
@test vel ≈ rvm[4:6]
vel = velocity(spk, 0, 1, [jd; jd], [0.0, 0.0])
@test vel ≈ mrvm[:,4:6]

st = state(spk, "mercury barycenter", jd)
@test st ≈ rvm
st = state(spk, 1, jd)
@test st ≈ rvm
st = state(spk, "ssb", "mercury barycenter", jd)
@test st ≈ rvm
st = state(spk, 0, "mercury barycenter", jd)
@test st ≈ rvm
st = state(spk, "ssb", 1, jd)
@test st ≈ rvm
st = state(spk, 0, 1, jd)
@test st ≈ rvm
st = state(spk, 0, 1, [jd; jd])
@test st ≈ mrvm
st = state(spk, "mercury barycenter", jd, 0.0)
@test st ≈ rvm
st = state(spk, 1, jd, 0.0)
@test st ≈ rvm
st = state(spk, "ssb", "mercury barycenter", jd, 0.0)
@test st ≈ rvm
st = state(spk, 0, "mercury barycenter", jd, 0.0)
@test st ≈ rvm
st = state(spk, "ssb", 1, jd, 0.0)
@test st ≈ rvm
st = state(spk, 0, 1, jd, 0.0)
@test st ≈ rvm
st = state(spk, 0, 1, [jd; jd], [0.0; 0.0])
@test st ≈ mrvm
