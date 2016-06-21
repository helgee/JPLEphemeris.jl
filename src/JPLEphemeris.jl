VERSION >= v"0.4" && __precompile__()

module JPLEphemeris

using Compat
import Compat: String

include("naif.jl")
include("daf.jl")
include("spk.jl")
include("ascii.jl")

end
