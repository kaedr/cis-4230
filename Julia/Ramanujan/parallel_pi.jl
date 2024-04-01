using Base.Threads
import("common_pi.jl")

function parallel_work()
    for k in threadid():nthreads():size
    end
end

function main()
end

main()
