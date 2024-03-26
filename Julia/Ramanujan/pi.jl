
struct Memo
    k::UInt;
    # (4n)! = (4k!) * 4n * 4(n-1) ... 4(k+1)
    four_k_factorial::BigInt
    # n!^4 = k!^4 * n^4 * (n-1)^4 ... (k+1)^4
    k_factorial_fourth::BigInt
    # 26390n = 26390k * 26390(n-k)
    k_26390::BigInt
    # 396^4n = 396^4k * 396^4(n-k)
    three_ninety_six_to_four_k::BigInt
end

# (4n)! = (4k!) * 4n * 4n-1 ... 4k+1
function next_four_k_factorial(memo::Memo, next_k::UInt)
    for i in (memo.k + 1):next_k
       fact_term = 4i
       memo.four_k_factorial *= fact_term
       memo.four_k_factorial *= fact_term - 1
       memo.four_k_factorial *= fact_term - 2
       memo.four_k_factorial *= fact_term - 3
    end
end

# n!^4 = k!^4 * n^4 * (n-1)^4 ... (k+1)^4
function next_k_factorial_fourth(memo::Memo, next_k::UInt)
    for i in (memo.k + 1):next_k
        memo.k_factorial_fourth *= i^4
    end
end

# 26390n = 26390k + 26390(n-k)
function next_k_26390(memo::Memo, next_k::UInt)
    memo.k_26390 += 26390 * next_k - memo.k
end

# 396^4n = 396^4k * 396^4(n-k)
function next_three_ninety_six_to_four_k(memo::Memo, next_k::UInt)

end

function main()
    precision = 0

    # Collect our target precision
    if length(ARGS) > 0
        precision = parse(Int, ARGS[1])
    end

    bit_precision = ceil(Int, precision * ( log( 10 ) / log( 2 )  ))
    setprecision(BigFloat, bit_precision)
    setrounding(BigFloat, RoundNearest)

    accumulator::BigFloat = 0.0
    k_memo = Memo(0, 1, 1, 0, 1)
end

main()
