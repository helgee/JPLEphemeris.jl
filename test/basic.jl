using Base.Test
using JPLEphemeris

# Reference value from CSPICE
rvm = [4.250906022073639e7,2.3501057648129586e7,8.158467467032234e6,-34.160008371029825,37.844059275357594,23.756128199757867]
mrvm = [rvm';rvm']

spk = SPK("$path/de430.bsp")
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
