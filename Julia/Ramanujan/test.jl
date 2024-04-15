include("common_pi.jl")

desired_precision::UInt = 1000000
bit_precision = ceil(Int, desired_precision * ( log( 10 ) / log( 2 )  ))
setprecision(BigFloat, bit_precision)
setrounding(BigFloat, RoundNearest)

function adjacent_k(sample_size::UInt)
    k_memo = Memo(0, 0, 1, 1, 1, 0.0)
    accumulator = BigFloat(0.0)

    for k in 0:sample_size
        advance_to_next_k(k_memo, k)
        calculate_term(k_memo)
        # Core.println("k: $k -> $(k_memo.term)")
        accumulator += k_memo.term
    end

end

function spaced_k(sample_size::UInt, spacing::UInt)
    k_memo = Memo(0, 0, 1, 1, 1, 0.0)
    accumulator = BigFloat(0.0)

    for k in 0:spacing:sample_size
        advance_to_next_k(k_memo, k)
        calculate_term(k_memo)
        # Core.println("k: $k -> $(k_memo.term)")
        accumulator += k_memo.term
    end

end

using Profile
using PProf
sample_size::UInt = 100

# println("Running adjacent_k")
# @time adjacent_k(sample_size)
# @time adjacent_k(sample_size)
# @profile adjacent_k(sample_size)
# Profile.print()
# Profile.clear()

# spacing::UInt = 4
# println("Running spaced_k 4")
# @time spaced_k(sample_size, spacing)
# @time spaced_k(sample_size, spacing)
# @profile spaced_k(sample_size, spacing)
# Profile.print()
# Profile.clear()

spacing::UInt = 24
println("Running spaced_k 24")
@time spaced_k(sample_size, spacing)
@time spaced_k(sample_size, spacing)
Profile.Allocs.@profile sample_rate=1 spaced_k(sample_size, spacing)
PProf.Allocs.pprof()
Profile.Allocs.clear()
