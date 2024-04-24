function threaded_sum(float_array)
    chunks = Iterators.partition(float_array, length(float_array) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn sum(chunk)
    end
    chunk_sums = fetch.(tasks)
    return sum(chunk_sums)
end


function main()
    array_length = 100000
    if length(ARGS) > 0
        array_length = parse(UInt, ARGS[1])
    end

    println()
    @time float_array = fill(1.0, array_length)

    println()
    println("Single threaded Sum:")
    @time sum(float_array)
    @time sum(float_array)
    @time sum(float_array)

    println()
    println("$(Threads.nthreads()) Threaded Sum:")
    @time threaded_sum(float_array)
    @time threaded_sum(float_array)
    @time threaded_sum(float_array)
end

main()
