function gauss_legendre(prec)
    setprecision(BigFloat, prec, base=10)
    GC.enable(false)
    println(precision(BigFloat, base=2))
    x::BigFloat = a::BigFloat = BigFloat(1, precision=precision(BigFloat, base=2))
    b::BigFloat = BigFloat(a / sqrt(BigFloat(2, precision=precision(BigFloat, base=2))))
    t::BigFloat = BigFloat(a / 4, precision=precision(BigFloat, base=2))
    y::BigFloat = a
    for _=0:ceil(BigInt, log2(prec))
        y, a = a::BigFloat, (a::BigFloat + b::BigFloat) / 2
        b = sqrt(b::BigFloat*y::BigFloat)
        t = t::BigFloat-(x::BigFloat * (y::BigFloat-a::BigFloat)^2)
        x = x::BigFloat*2
    end
    (a + b)^2 / (4 * t)
end

@time begin
    my_pi = gauss_legendre(1000000)
end

println(Float64((my_pi::BigFloat-BigFloat(π))::BigFloat*1000000))
