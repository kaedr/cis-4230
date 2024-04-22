@everywhere using Base.Threads
@everywhere using Distributed

@everywhere function greet(greeted, greeter)
    println("Hello $greeted, I'm $greeter")
end

hello(name) = greet(name, "Julia")

function threaded_hello(name)
    @threads for _ in 1:nthreads()
        greet(name, "Thread $(threadid())")
    end
end

@everywhere function distributed_hello(name)
    @threads for _ in 1:nthreads()
        greet(name, "Process $(myid()), Thread $(threadid())")
    end
end

if !isinteractive()
    name = "User"
    if length(ARGS) > 0
        name = ARGS[1]
    end

    if nprocs() > 1
        @everywhere name = $name
        return @everywhere distributed_hello(name)
    end

    nthreads() > 1 && return threaded_hello(name)

    hello(name)
end
