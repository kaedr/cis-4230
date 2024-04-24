using Printf

# For correct initial values, construct as:
# Memo(0, 0, 1, 1, 1, 0.0)
mutable struct Memo
    k::UInt;
    # 26390n = 26390k * 26390(n-k)
    k_26390::BigInt
    # 396^4n = 396^4k * 396^4(n-k)
    three_ninety_six_to_four_k::BigInt
    # (4n)! = (4k!) * 4n * 4(n-1) ... 4(k+1)
    four_k_factorial::BigInt
    # n!^4 = k!^4 * n^4 * (n-1)^4 ... (k+1)^4
    k_factorial_fourth::BigInt
    # Avoid reinitializing the term each loop
    term::BigFloat
end

default_memo() = Memo(0, 0, 1, 1, 1, 0.0)

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
        memo.k_factorial_fourth *= BigInt(i)^4
    end
end

# 26390n = 26390k + 26390(n-k)
function next_k_26390(memo::Memo, next_k::UInt)
    memo.k_26390 += 26390 * (next_k - memo.k)
end

# 396^4n = 396^4k * 396^4(n-k)
function next_three_ninety_six_to_four_k(memo::Memo, next_k::UInt)
    # In initially forgot to make this a bigint, which worked fine for Serial
    # But broke as soon as the term overflowed an integer
    memo.three_ninety_six_to_four_k *= BigInt(396) ^ 4(next_k - memo.k)
end

function advance_to_next_k(memo::Memo, next_k::UInt)
    next_four_k_factorial(memo, next_k)
    next_k_factorial_fourth(memo, next_k)
    next_k_26390(memo, next_k)
    next_three_ninety_six_to_four_k(memo, next_k)
    memo.k = next_k
    # println(memo)
end

function calculate_term(memo::Memo)
    # Set to 26390k + 1103
    memo.term = memo.k_26390
    memo.term += 1103
    # Divide by 396^4k
    memo.term /= memo.three_ninety_six_to_four_k
    # Multiply by 4k!
    memo.term *= memo.four_k_factorial
    # Divide by k!^4
    memo.term /= memo.k_factorial_fourth
end

function check_results(accumulator::BigFloat, target_precision::UInt)
    # println("Checking: $accumulator")
    # Multiplying by 10 ^ target_precision shifts our difference a number of positions left
    # equal to our desired_precision, making it easier to see the scope of any discrepency
    # The downcast to Float64 is to make it display nicely
    difference = Float64((BigFloat(π) - accumulator)::BigFloat * BigInt(10) ^ target_precision)
    println("Difference at $target_precision digits is $difference")
end

function check_string_results(accumulator::BigFloat, target_precision::UInt)
    pi_file_string = read("pi.txt", String)
    pi_calc_string = string(accumulator)

    accuracy = 0
    for (index, (ref, calced)) in enumerate(zip(pi_file_string, pi_calc_string))
        # We set this before checking accuracy because the nature of rounding means
        # We will be wrong about 50% of the time on the last digit
        accuracy = index
        if ref == calced
            if index >= target_precision
                break
            end
        else
            break
        end
    end
    println("Calculation accurate to $(accuracy) digits")
end
