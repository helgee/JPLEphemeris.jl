JPLEphemeris.jl
===============

[![Travis Status](https://travis-ci.org/helgee/JPLEphemeris.jl.png)](https://travis-ci.org/helgee/JPLEphemeris.jl)
[![AppVeyor status](https://ci.appveyor.com/api/projects/status/7pt2vy8wulix06jk?svg=true)](https://ci.appveyor.com/project/helgee/jplephemeris-jl)
[![PkgEval v4](http://pkg.julialang.org/badges/JPLEphemeris_0.4.svg)](http://pkg.julialang.org/?pkg=JPLEphemeris)
[![PkgEval v5](http://pkg.julialang.org/badges/JPLEphemeris_0.5.svg)](http://pkg.julialang.org/?pkg=JPLEphemeris)

The [JPL Development Ephemerides][jpl] are the results of simulations of the Solar System used for spacecraft navigation and astronomical purposes. They are published as [SPK kernel files][spk] which contain sets of Chebyshev polynomial coefficients with which the position and velocity of the Solar System's planets can be interpolated with high precision for all dates covered by the ephemeris.

This package provides functionality to read SPK files and compute the position and velocity of the planets and minor bodies directly from Julia.

## Installation

The package can be installed through Julia's package manager.

```julia
Pkg.add("JPLEphemeris")
```

## Usage

```julia
using JPLEphemeris

# Load the DE430 SPK kernel
spk = SPK("de430.bsp")

# 2016-01-01T00:00 in Julian days
jd = Dates.datetime2julian(DateTime(2016,1,1,0,0,0))

# Position of Earth's barycenter w.r.t. the Solar System's barycenter at 2016-01-01T00:00
# [km]
pos = position(spk, "earth barycenter", jd)

# Velocity of Earth w.r.t. Earth's barycentre at 2016-01-01T00:00
# [km/s]
vel = velocity(spk, "earth barycenter", "earth", jd)

# Compute the state vector (position and velocity) of Earth's barycenter (NAIF ID: 3)
# w.r.t. to the Solar System's barycenter (NAIF ID: 0) for a range of Julian days
st = state(spk, 0, 3, jd:jd+100)

# Two-part Julian dates (day number and fraction) can be used for higher precision.
# For example for 2016-01-01T12:00:
st = state(spk, 0, 3, jd, 0.5)
```

## ASCII API

This package also provides an older API for working with ephemeris data from ASCII tables which are converted to [JLD files][jld].
However using the SPK API is recommended.

### Usage

```julia
using JPLEphemeris.ASCII

# Download and convert ephemeris DE421
getephem(421)

# Load ephemeris DE421
eph = Ephemeris(421)

# 2014-01-01T00:00 in Julian days
jd = 2456658.5

# Position of Mercury w.r.t. the Solar System's barycentre at 2014-01-01T00:00
# [km]
pos = position(eph, "mercury", jd)

# Velocity of Mercury w.r.t. the Solar System's barycentre at 2014-01-01T00:00
# [km/day]
vel = velocity(eph, "mercury", jd)

# Complete state vector (position and velocity) for a range of Julian days
st = state(eph, "mercury", jd:jd+100)

# The ephemeris also contains the set of constants with which it was calculated
# e.g. the Astronomical Unit (AU)
au = eph.constants["AU"]

# Close ephemeris file
close(eph)
```

### Available Celestial Bodies

* Mercury: `"mercury"`
* Venus: `"venus"`
* Earth-Moon system's barycenter: `"earthmoon"`
* Mars: `"mars"`
* Jupiter: `"jupiter"`
* Saturn: `"saturn"`
* Uranus: `"uranus"`
* Neptune: `"neptune"`
* Pluto: `"pluto"`
* Moon (geocentric): `"moon"`
* Sun: `"sun"`
* Librations: `"librations"`
* Nutations: `"nutations"`

## Validation

The package's test suite uses test files provided by JPL to confirm that the results are correct.

## Acknowlegements
Most of this is based on the excellent [jplehem][jplephem] library by [Brandon Rhodes][br].
Please use it if you need similar functionality in Python.

[jpl]: http://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
[jplephem]: https://github.com/brandon-rhodes/python-jplephem
[jld]: https://github.com/JuliaLang/JLD.jl
[br]: https://github.com/brandon-rhodes
[spk]: http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/
