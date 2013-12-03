JPLEphemeris.jl
===============

This package provides a convenient way to use the [JPL Development Ephemerides][jpl] from Julia. [JLD][hdf5]

## Installation

```julia
Pkg.add("JPLEphemeris")
```

## Usage

```julia
using JPLEphemeris

eph = Ephemeris()

# 2014-01-01T00:00 in Julian days
jd = 2456658.5

# Position of Mercury w.r.t. the Solar System's barycentre
pos = position(eph, "mercury", jd)

# Velocity of Mercury w.r.t. the Solar System's barycentre
vel = velocity(eph, "mercury", jd)

# Complete state vector (position and velocity)
st = state(eph, "mercury", jd)
```

## Acknowlegements
Most of this is based on the excellent [jplehem][jplephem] library by Brandon Rhodes.
Please use it if you need similar functionality in Python.

[jpl]: http://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
[jplephem]: https://github.com/brandon-rhodes/python-jplephem
[hdf5]: https://github.com/timholy/HDF5.jl
