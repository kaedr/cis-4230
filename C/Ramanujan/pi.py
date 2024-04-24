import math

#computes pi value by Ramanujan formula
#source : https://www.maa.org/sites/default/files/pdf/upload_library/22/Chauvenet/BorweinBorweinBailey.pdf
def pi_formula():
    sum = 0
    n = 0
    i = (math.sqrt(8))/9801

    while True:
        four_n_factorial = math.factorial(4*n)
        n_factorial_fourth = pow(math.factorial(n),4)
        n_26390_plus_1103 = (26390*n+1103)
        pow_term = pow(396,4*n)
        print(f"{n_26390_plus_1103} / {pow_term} * {four_n_factorial} / {n_factorial_fourth}")
        tmp = n_26390_plus_1103
        print(f"term value is: {tmp:.30g}")
        tmp /= pow_term
        print(f"term value is: {tmp:.30g}")
        tmp *= four_n_factorial
        print(f"term value is: {tmp:.30g}")
        tmp /= n_factorial_fourth
        print(f"term value is: {tmp:.30g}")
        sum +=tmp
        print(f"accumulator value is: {sum:.30g}")

        if(abs(tmp) < 1e-30):
            break
        n += 1

    sum *= i

    return(1/sum)

pi = pi_formula()
print(f"Pi value using Ramanujan Formula : {pi:.30g}")
