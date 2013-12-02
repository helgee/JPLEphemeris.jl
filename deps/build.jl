using JPLEphemeris.Parser

const DENUM = "421"
const FTP = "ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de$DENUM"
const HEADER = "header.421"
const DATAFILES = ["ascp1900.421", "ascp2050.421"]
const TESTFILE = "testpo.421"
const OUTFILE = "de421.jld"

if ~isfile(HEADER)
    download("$FTP/$HEADER", "$HEADER")
end
if ~isfile(TESTFILE)
    download("$FTP/$TESTFILE", "$TESTFILE")
end
for f in DATAFILES
    if ~isfile(f)
        download("$FTP/$f", "$f")
    end
end
if ~isfile(OUTFILE)
    readascii(HEADER, DATAFILES, OUTFILE)
end

# Remove ASCII files
rm(HEADER)
for f in DATAFILES
    rm(f)
end
