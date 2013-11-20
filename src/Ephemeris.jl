module Ephemeris

function poly(a:Vector, x:Float64)
    y = 0.0
    for i = 1:length(a)
        y += a[i]*x^i
    end
    return y
end

end
