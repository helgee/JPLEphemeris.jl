module Ephemerides

type Ephemeris
    startepoch::Float64
    finalepoch::Float64
end

function readascii(header)
    startepoch = 0.0
    finalepoch = 0.0
    firstline = 0
    lastline = 0
    l = open(readlines, header)
    for i = 1:length(l)
        if beginswith(l[i], "GROUP   1030")
            startepoch, finalepoch, _ = split(l[i+2])
        elseif beginswith(l[i], "GROUP   1041")
            n = int(l[i+2])
            firstline = i+2
            lastline = i+2+div(n,3)
        end
    end
    b = zeros(Float64, lastline-firstline, 3)
    for i = 1:size(b)[1]
        println(l[i+firstline])
        b[i,:] = float(split(replace(l[i+firstline],"D","e")))
    end
    return startepoch, finalepoch, b
end

function poly(a::Vector, x::Float64)
    y = 0.0
    for i = 1:length(a)
        y += a[i]*x^i
    end
    return y
end

end
