using Distributed
using Base.Threads

const HOSTS = ["node1","node2","node3","node4"]
const THREADS_PER_HOST::UInt = nthreads()

# Set up remote processors
# Note that you need to run addprocs before any
# @everywhere declarations in order for those items to
# be pushed out to the added processors
addprocs(HOSTS, exeflags=["-t $THREADS_PER_HOST"])
@everywhere using Base.Threads
@everywhere include("common_pi.jl")
@everywhere include("distributed_funcs.jl")

function main()
    desired_precision::UInt = 100

    # Collect our target precision
    if length(ARGS) > 0
        desired_precision = parse(UInt, ARGS[1])
    end

    bit_precision = ceil(Int, desired_precision * ( log( 10 ) / log( 2 )  ))
    # We need precision set everywhere
    @everywhere bit_precision = $bit_precision
    @everywhere setprecision(BigFloat, bit_precision)
    @everywhere setrounding(BigFloat, RoundNearest)

    factor::BigFloat = sqrt(big(8)) / 9801


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

    # println(accumulator)

    check_results(accumulator, desired_precision)
    check_string_results(accumulator, desired_precision)
end

main()
