ftp = "ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii"
denums = ["431","430","424","423","422","421","418","414","413","410","406","405","403","202","200","102"]
headers = Dict(denums[3:end],["header.$d" for d in denums[3:end]])
merge!(headers, ["431"=>"header.431_572","430"=>"header.430_572"])

datafiles = ["430"=>["ascp$i.430" for i = 1550:100:2550],
"431"=>[[@sprintf("ascm%05d.431", i) for i = 1000:1000:13000], [@sprintf("ascp%05d.431", i) for i = 0:1000:16000]],
"424"=>[[@sprintf("ascm%04d.424", i) for i = 100:100:3000], [@sprintf("ascp%04d.424", i) for i = 0:100:2900]],
"423"=>["ascp$i.423" for i = 1800:50:2150],
"422"=>[[@sprintf("ascm%04d.422", i) for i = 100:100:3000], [@sprintf("ascp%04d.422", i) for i = 0:100:2900]],
"421"=>["ascp1900.421", "ascp2050.421"],
"418"=>["ascp1900.418"],
"414"=>["ascp$i.414" for i = 1600:100:2100],
"413"=>["ascp$i.413" for i = 1900:25:2025],
"410"=>["ascp$i.410" for i = 1960:20:2000],
"406"=>[[@sprintf("ascm%04d.406", i) for i = 100:100:3000], [@sprintf("ascp%04d.406", i) for i = 0:100:2900]],
"405"=>["ascp$i.405" for i = 1600:20:2200],
"403"=>["ascp$i.403" for i = 1600:100:2100],
"202"=>["ascp$i.202" for i = 1900:50:2150],
"200"=>["ascp$i.200" for i = 1600:20:2160],
"102"=>[[@sprintf("ascm%04d.102", i) for i = 200:300:1400], [@sprintf("ascp%04d.102", i) for i = 100:300:2800]]]

function download(ftp, file, out)
    errorstring = "\nPlease download all files from '$ftp' manually, place them in '$PATH' and re-run 'getephem'.\n"
    try
        if OS_NAME == :Windows
            run(`where curl`)
        else
            run(`which curl`)
        end
    catch
        error("Curl could not be found. $errorstring")
    end

    try
        run(`curl $ftp/$file -o $out`)
    catch
        error("Curl could not connect to the server. $errorstring")
    end
end

function getephem(denum; force=false)
    denum = string(denum)
    if ~haskey(datafiles, denum)
        error("Unknown ephemeris 'DE$denum'.")
    end
    outfile = "$PATH/de$denum.jld"
    if isfile(outfile) && force
        rm(outfile)
    end
    if ~isfile(outfile)
        testfile = "testpo.$denum"
        testlocal = "$PATH/$testfile"
        header = headers[denum]
        headerlocal = "$PATH/$header"
        data = datafiles[denum]
        datalocal = ["$PATH/$d" for d in data] 
        println("Building ephemeris DE$denum.")
        if ~isfile("$PATH/$header") || force
            download("$ftp/de$denum", header, headerlocal)
        end
        if ~isfile("$PATH/$testfile")
            download("$ftp/de$denum", testfile, testlocal)
        end
        for (f,d) in zip(data, datalocal)
            if ~isfile(d) || force
                download("$ftp/de$denum", f, d)
            end
        end
        readascii(headerlocal, datalocal, outfile)
        println("Removing ASCII files.")
        rm(headerlocal)
        for f in datalocal
            rm(f)
        end
    end
end

function rmephem(denum)
    rm("$PATH/de$denum.jld")
    rm("$PATH/testpo.$denum")
end

function readascii(header, datafiles, outfile)
    startepoch = 0.0
    finalepoch = 0.0
    ncoeff = 0
    constantnames = String[]
    constantvalues = Float64[]
    ind = zeros(Int64,3,13)
    intervall = zeros(13)
    dtable = Array(Vector{Float64},13)
    for i = 1:13
        dtable[i] = Float64[]
    end

    println("Parsing header file.")
    l = open(readlines, header)
    ncoeff = int(split(l[1])[4])
    for i = 1:length(l)
        if beginswith(l[i], "GROUP   1030")
            startepoch, finalepoch = float(split(l[i+2])[1:2])
        elseif beginswith(l[i], "GROUP   1040")
            n = int(l[i+2])
            firstline = i+3
            lastline = i+3+div(n,10)
            for j = firstline:lastline
                append!(constantnames, split(l[j]))
            end
        elseif beginswith(l[i], "GROUP   1041")
            n = int(l[i+2])
            firstline = i+3
            lastline = i+3+div(n,3)
            for j = firstline:lastline
                append!(constantvalues, float(split(replace(l[j],"D","e"))))
            end
        elseif beginswith(l[i], "GROUP   1050")
            firstline = i+2
            lastline = i+4
            for j = firstline:lastline
                ind[j-firstline+1,:] = int(split(l[j]))
            end
        end
    end
    constants = Dict(constantnames, constantvalues)

    file = jldopen(outfile, "w")
    try
        @write file constants
        @write file startepoch
        @write file finalepoch
        @write file ind
        for i = 1:length(datafiles)
            println("Parsing data file $i/$(length(datafiles)).")
            parsedatafile(datafiles[i], ncoeff, ind, dtable, file)
        end
        for i = 1:length(dtable)
            sort!(dtable[i])
            intervall[i] = dtable[i][2] - dtable[i][1]
        end
        @write file dtable
        @write file intervall
    catch
        close(file)
        rm(outfile)
        rethrow()
    finally
        close(file)
    end
end

function parsedatafile(datafile, ncoeff, ind, dtable, file)
    intervall = r"^\s+[0-9]+\s+[0-9]+"
    coeff = Float64[]
    intervalls = zeros(size(ind)[2])
    l = open(readlines, datafile)
    counter = 0
    for i = 1:length(l)
        if ismatch(intervall, l[i])
            counter += 1
            date, finaldate = float(split(replace(l[i+1],"D","e"))[1:2])
            firstline = i+1
            lastline = i+iceil(ncoeff/3)
            for j = firstline:lastline
                append!(coeff, float(split(replace(l[j],"D","e"))))
            end
            for j = 1:13
                dt = (finaldate - date)/ind[3,j]
                if ~in(date, dtable[j]) 
                    push!(dtable[j], date)
                elseif in(date, dtable[j]) 
                    continue
                end
                if j == 12
                    n = 2
                else
                    n = 3
                end
                offset = ind[2,j]*n
                i1 = ind[1,j]
                i2 = i1 + offset - 1
                write(file, "$j-$date", reshape(coeff[i1:i2], ind[2,j], n))
                if ind[3,j] != 1
                    dt = (finaldate - date)/ind[3,j]
                    for k = 1:ind[3,j]-1
                        date1 = date + dt*k
                        push!(dtable[j], date1)
                        i1k = i1 + k*offset
                        i2k = i1 + (k+1)*offset - 1
                        write(file, "$j-$date1",
                            reshape(coeff[i1k:i2k], ind[2,j], n))
                    end
                end
            end
            empty!(coeff)
        end
    end
end
