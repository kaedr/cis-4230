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
        print(f"{four_n_factorial} / {n_factorial_fourth} * {n_26390_plus_1103} / {pow_term}")
        tmp = i*(four_n_factorial/n_factorial_fourth*n_26390_plus_1103/pow_term)
        print(f"iteration value is: {tmp}")
        sum +=tmp

        if(abs(tmp) < 1e-15):
            break
        n += 1


    return(1/sum)

print("Pi value using Ramanujan Formula : ",pi_formula())
