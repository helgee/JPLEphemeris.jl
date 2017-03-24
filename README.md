JPLEphemeris.jl
===============

[![Travis Status](https://travis-ci.org/JuliaAstrodynamics/JPLEphemeris.jl.png)](https://travis-ci.org/JuliaAstrodynamics/JPLEphemeris.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/2d9y5mcek1qxggg0?svg=true)](https://ci.appveyor.com/project/JuliaAstrodynamics/jplephemeris-jl)
[![PkgEval v5](http://pkg.julialang.org/badges/JPLEphemeris_0.5.svg)](http://pkg.julialang.org/?pkg=JPLEphemeris)

The [JPL Development Ephemerides][jpl] are the results of simulations of the Solar System used for spacecraft navigation and astronomical purposes. They are published as [SPK kernel files][spk] which contain sets of Chebyshev polynomial coefficients with which the position and velocity of the Solar System's planets can be interpolated with high precision for all dates covered by the ephemeris.

This package provides functionality to read SPK files and compute the position and velocity of the planets directly from Julia.

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

# List the available segments
print_segments(spk)

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

The ASCII API that was originally included with this package has been moved to [LegacyEphemeris.jl][legacy].

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
[legacy]: https://github.com/helgee/LegacyEphemeris.jl
