using JPLEphemeris
using Base.Test
import Base.source_path

path = dirname((source_path()))

if ~isfile("$path/header.421")
    run(`curl -o $path/header.421 ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de421/header.421`)
end

if ~isfile("$path/ascp1900.421")
    run(`curl -o $path/ascp1900.421 ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de421/ascp1900.421`)
end

if ~isfile("$path/ascp2050.421")
    run(`curl -o $path/ascp2050.421 ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de421/ascp2050.421`)
end

header = "$path/header.421"
datafiles = ["$path/ascp1900.421", "$path/ascp2050.421"]
outfile = "$path/de421.jld"

coeff = readascii(header, datafiles, outfile)
