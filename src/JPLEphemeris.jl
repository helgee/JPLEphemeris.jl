VERSION >= v"0.4" && __precompile__()

module JPLEphemeris

include("naif.jl")
include("daf.jl")
include("spk.jl")
include("ascii.jl")

end
