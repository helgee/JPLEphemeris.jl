JPLEphemeris.jl
===============

[![Build Status](https://travis-ci.org/helgee/JPLEphemeris.jl.png)](https://travis-ci.org/helgee/JPLEphemeris.jl)
[![Coverage Status](https://img.shields.io/coveralls/helgee/JPLEphemeris.jl.svg)](https://coveralls.io/r/helgee/JPLEphemeris.jl)

The [JPL Development Ephemerides][jpl] are the results of simulations of the Solar System used for spacecraft navigation and astronomical purposes. They are published as sets of Chebyshev polynomial coefficients with which the position and velocity of the Solar System's planets can be interpolated with high precision for all dates covered by the ephemeris.
This package provides functionality to convert the data from the original ASCII tables to the HDF5-based [JLD format][jld] for efficient access and routines for the computation of the planet's state vectors directly from Julia.

## Requirements

CURL must be installed for the automatic download of ephemeris files.

## Installation

The package can be installed through Julia's package manager like shown below:

```julia
Pkg.add("JPLEphemeris")
```
A standard ephemeris (currently DE430) will be downloaded and compiled automatically during the installation. Additional ephemerides can be added later.

## Usage

```julia
using JPLEphemeris

# Initialize the standard ephemeris DE430
eph = Ephemeris()

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

The `position`, `velocity` and `state` functions can also be called with arrays or ranges of Julian days.

Another example is shown in this IJulia [notebook][notebook].

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

### Adding Additional Ephemerides

All ephemerides available at ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/ can be downloaded and converted automatically with the `getephem` command.

```julia
using JPLEphemeris

# e.g. DE421
getephem(421)
eph = Ephemeris(421)
```

If CURL is not available the user will be prompted to fetch the files manually and to re-run `getephem` subsequently.
Should the build process fail, e.g. due to corrupted ASCII files, re-downloading of the data files can be forced via the `force` keyword parameter.

```julia
getephem(421, force=true)
```

Ephemeris files can be removed with `rmephem(421)`.

## Validation

The package's test suite uses test files provided by JPL to confirm that the results are correct for all installed ephemerides.

## Acknowlegements
Most of this is based on the excellent [jplehem][jplephem] library by [Brandon Rhodes][br].
Please use it if you need similar functionality in Python.

[jpl]: http://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
[jplephem]: https://github.com/brandon-rhodes/python-jplephem
[jld]: https://github.com/timholy/HDF5.jl
[notebook]: http://nbviewer.ipython.org/github/helgee/JPLEphemeris.jl/blob/master/JPLEphemeris-Earth_Mars-2014.ipynb
[br]: https://github.com/brandon-rhodes
