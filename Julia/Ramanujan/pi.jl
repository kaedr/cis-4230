include("common_pi.jl")


function main()
    desired_precision::UInt = 100

    # Collect our target precision
    if length(ARGS) > 0
        desired_precision = parse(UInt, ARGS[1])
    end

    bit_precision = ceil(Int, desired_precision * ( log( 10 ) / log( 2 )  ))
    setprecision(BigFloat, bit_precision)
    setrounding(BigFloat, RoundNearest)

    accumulator::BigFloat = 0.0
    k_memo = Memo(0, 0, 1, 1, 1, big"0.0")
    factor::BigFloat = sqrt(big(8)) / 9801
    # I found that each iteration adds closer to 7 digit of precision
    iterations::UInt = ceil(desired_precision / 7)
    println("Making $iterations iterations")

    start_time = time()
    for k in 0:iterations
        advance_to_next_k(k_memo, k)
        calculate_term(k_memo)
        accumulator += k_memo.term
        # println(k_memo)
        # println(accumulator)
        if k > 0 && k % 2000 == 0
            checkpoint = time() - start_time
            println("$k iterations complete in $checkpoint seconds")
        end
    end
    checkpoint = time() - start_time
    println("$iterations iterations complete in $checkpoint seconds")

    # multiply the factor into our accumulator
    accumulator *= factor
    # account for the 1/pi thing
    accumulator = 1 / accumulator

    check_results(accumulator, desired_precision)
    check_string_results(accumulator, desired_precision)
end

main()
