VERSION >= v"0.4" && __precompile__()

module JPLEphemeris

using Compat
import Compat: ASCIIString

include("naif.jl")
include("daf.jl")
include("spk.jl")

end
