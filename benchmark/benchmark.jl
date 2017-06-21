using AstroDynBase
using BenchmarkTools
using Distributions
using JPLEphemeris

srand(430)

function runner(times)
    arr = zeros(6)
    state!(arr, spk, seg, 1.0, times[1])
    total = 0.0
    for t in times
        time = @elapsed state!(arr, spk, seg, 1.0, t)
        total += time
    end
    print_time(total / n)
end

function print_time(time)
    if time < 1e-6
        println("Time: $(time*1e9) ns")
    elseif time < 1e-3
        println("Time: $(time*1e6) μs")
    elseif time < 1.0
        println("Time: $(time*1e3) ms")
    else
        println("Time: $time s")
    end
end

spk = SPK(joinpath(Pkg.dir("JPLEphemeris"), "deps", "de430.bsp"))
seg = spk.segments[0][3]
n = 10_000_000

first = seg.firstdate
last = seg.lastdate
linear = collect(linspace(first, last, n))
d = Truncated(Normal(), first, last)
random = rand(d, n)

runner(linear)
runner(random)
# @profiler state.(spk, seg, random)
# Profile.clear()
# @profile state.(spk, seg, random)
# Profile.print()
