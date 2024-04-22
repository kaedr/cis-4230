# Crash Course in Julia

## Introduce the REPL

## Introductory information
1. [Julia Docs](https://docs.julialang.org/en/v1/)
2. [Arbitrary Precision Arithmetic](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic)
    1. [BigInt](https://github.com/JuliaLang/julia/blob/b5bfd83a3d0ee55f27fbf18dfb9761a3f284fd99/base/gmp.jl)
    2. [BigFloat](https://github.com/JuliaLang/julia/blob/b5bfd83a3d0ee55f27fbf18dfb9761a3f284fd99/base/mpfr.jl)
3. [Mathmatical Operations](https://docs.julialang.org/en/v1/manual/mathematical-operations/)
4. Symbology (√ and π)
    1. [Whence Commeth pi](https://www.mpfr.org/algorithms.pdf)
    2. [Brent-Salamin formula](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_algorithm)


## Hello

1. basic hello `julia hello.jl`
2. [threaded](https://docs.julialang.org/en/v1/manual/multi-threading/) hello
    1. static number of threads `julia -t 2 hello.jl`
    2. dynamic number of threads `julia -t auto hello.jl`
2. [distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) hello
    1. local processes `julia -p 2 -t 2 hello.jl`
    2. remote processes `julia --machine-file=hosts.txt -t 4 hello.jl`

## Dig into pi in Julia

## Performance analysis
- `@profile` macro
    ```julia
    using Profile
    @profile some_function()
    Profile.print()

    # If you want to get rid of accumulated results
    Profile.clear()
    ```
- `@time` macro
- `Profile.Allocs.@profile` + PProf.jl
    ```Julia
    Profile.Allocs.@profile sample_rate=1 some_function()
    PProf.Allocs.pprof()
    Profile.Allocs.clear()
    ```

Julia [tips](https://docs.julialang.org/en/v1/manual/performance-tips/)

