using Base.Threads

function thread_work(start_point::UInt, interval::UInt, iterations::UInt)
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

function process_work(num_procs::UInt, num_threads::UInt, iterations::UInt)
    # Account for the first worker being ID 2
    start_point = (Distributed.myid() - 2) * num_threads
    interval = num_procs * num_threads
    accumulator::BigFloat = 0.0
    workers = Array{Task}(undef, num_threads)

    # Account for Julia's inclusive ranges
    for i::UInt in 1:num_threads
        # start_point is the beginning of this processors work block
        # the -1 is to offset for Julia's 1 based indexing
        thread_start_point = start_point + i - 1
        # println("Thread $i starting at $thread_start_point")
        workers[i] = Threads.@spawn thread_work(thread_start_point, interval, iterations)
    end

    for i in 1:num_threads
        accumulator += fetch(workers[i])
    end

    return accumulator
end

