using Distributed
@everywhere include("common_pi.jl")

@everywhere function thread_work(tid::UInt, iterations::UInt)
    k_memo = Memo(0, 0, 1, 1, 1, 0.0)
    accumulator = BigFloat(0.0)
    start_time = time()
    # Have to do tid - 1 because of Julia's 1 based everything
    for k in tid-1:nthreads():iterations
        advance_to_next_k(k_memo, k)
        calculate_term(k_memo)
        # Core.println("k: $k -> $(k_memo.term)")
        accumulator += k_memo.term
        # This will only be run by one thread because only one will run the item
        # that is evenly divisble by 2000
        if k > 0 && k % 2000 == 0
            checkpoint = time() - start_time
            Core.println("$k iterations complete in $checkpoint seconds")
        end
    end
    # Thread id 1 gives the final report
    if tid == 1
        checkpoint = time() - start_time
        Core.println("$iterations iterations complete in $checkpoint seconds")
    end
    return accumulator
end

function main()
    desired_precision::UInt = 100

    # Collect our target precision
    if length(ARGS) > 0
        desired_precision = parse(UInt, ARGS[1])
    end

    bit_precision = ceil(Int, desired_precision * ( log( 10 ) / log( 2 )  ))
    setprecision(BigFloat, bit_precision)
    setrounding(BigFloat, RoundNearest)

    factor::BigFloat = sqrt(big(8)) / 9801
    # I found that each iteration adds closer to 7 digit of precision
    iterations::UInt = ceil(desired_precision / 7)
    nodes = nprocs()
    println("Making $iterations iterations across $nodes processes")

    accumulator::BigFloat = 0.0
    workers = Array{Future}(undef, num_procs)

    start_time = time()

    # Do work here

    checkpoint = time() - start_time
    println("$iterations iterations complete in $checkpoint seconds")

    # multiply the factor into our accumulator
    accumulator *= factor
    # account for the 1/pi thing
    accumulator = 1 / accumulator

    # println(accumulator)

    check_results(accumulator, desired_precision)
end

main()
