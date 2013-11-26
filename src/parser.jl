module Parser

using HDF5, JLD

export readascii

function readascii{T<:String}(header::String, datafiles::Vector{T}, outfile::String)
    startepoch = 0.0
    finalepoch = 0.0
    ncoeff = 0
    constantnames = String[]
    constantvalues = Float64[]
    ind = zeros(Int64,3,13)
    intervall = zeros(Float64, 13)
    dtable = Array(Vector{Float64},13)
    for i = 1:13
        dtable[i] = Float64[]
    end

    # Parse header file
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

    # Save data to JLD file
    file = jldopen(outfile, "w")
    try
        @write file constants
        @write file startepoch
        @write file finalepoch
        @write file ind
        # Parse data files
        for i = 1:length(datafiles)
            parsedatafile(datafiles[i], ncoeff, ind, dtable, file)
        end
        for i = 1:length(dtable)
            sort!(dtable[i])
            intervall[i] = dtable[i][2] - dtable[i][1]
        end
        @write file dtable
        @write file intervall
    finally
        close(file)
    end
end

function parsedatafile(datafile, ncoeff, ind, dtable, file)
    intervall = r"^\s+[0-9]+\s+[0-9]+"
    coeff = Float64[]
    intervalls = zeros(Float64, size(ind)[2])
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

end
