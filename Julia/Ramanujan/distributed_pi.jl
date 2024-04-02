using Distributed
@everywhere using Base.Threads
@everywhere include("common_pi.jl")

const HOSTS = ["node1","node2","node3","node4"]
const THREADS_PER_HOST::UInt = nthreads()

@everywhere function thread_work(start_point::UInt, interval::UInt, iterations::UInt)
    k_memo = Memo(0, 0, 1, 1, 1, 0.0)
    accumulator = BigFloat(0.0)
    start_time = time()

    for k in start_point:interval:iterations
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

    return accumulator
end

@everywhere function process_work(num_procs::UInt, num_threads::UInt, iterations::UInt)
    # Account for Julia's 1 based indexing
    start_point = Distributed.myid() - 1
    interval = num_procs * num_threads
    accumulator::BigFloat = 0.0
    workers = Array{Task}(undef, num_threads)

    # Account for Julia's inclusive ranges
    for i::UInt in start_point:start_point + num_threads -1
        workers[i] = Threads.@spawn thread_work(i, interval, iterations)
    end

    for i in 1:num_threads
        accumulator += fetch(workers[i])
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

    # Set up remote processors
    # addprocs(HOSTS, exeflags=["-t $num_threads"])

    # I found that each iteration adds closer to 7 digit of precision
    iterations::UInt = ceil(desired_precision / 7)
    nodes::UInt = nworkers()
    println("Making $iterations iterations across $nodes processors with $THREADS_PER_HOST threads each")

    accumulator::BigFloat = 0.0
    result_futures = Array{Future}(undef, nprocs())

    start_time = time()

    # Distribute work to processors (local and remote)
    for worker in workers()
        result_futures[worker] = Distributed.@spawnat worker process_work(nodes, THREADS_PER_HOST, iterations)
    end

    for worker in workers()
        accumulator += fetch(result_futures[worker])
    end

    checkpoint = time() - start_time
    println("$iterations iterations complete in $checkpoint seconds")

    # multiply the factor into our accumulator
    accumulator *= factor
    # account for the 1/pi thing
    accumulator = 1 / accumulator

    println(accumulator)

    check_results(accumulator, desired_precision)
end

main()
