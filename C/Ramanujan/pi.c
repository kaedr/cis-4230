
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gmp.h>
#include <mpfr.h>

#define DECIMAL_PRECISION    10000   // Digits of precision

// A struct for memoization of previous calculations
struct Memo {
    // Below are expressed the calculations of these terms
    // for any future k (which we'll call n)
    unsigned long k;
    // (4n)! = (4k!) * 4n * 4(n-1) ... 4(k+1)
    mpz_t four_k_factorial;
    // n!^4 = k!^4 * n^4 * (n-1)^4 ... (k+1)^4
    mpz_t k_factorial_fourth;
    // 26390n = 26390k * 26390(n-k)
    mpz_t k_26390;
    // 396^4n = 396^4k * 396^4(n-k)
    mpz_t three_ninety_six_to_four_k;
};

// Intialize the Memo to a useful starting state
void memo_init(struct Memo *k_memo) {
    k_memo->k = 0;
    mpz_init_set_ui( k_memo->four_k_factorial, 1UL );
    mpz_init_set_ui( k_memo->k_factorial_fourth, 1UL );
    mpz_init_set_ui( k_memo->k_26390, 0UL );
    mpz_init_set_ui( k_memo->three_ninety_six_to_four_k, 1UL );
}

void memo_destroy(struct Memo *k_memo) {
    mpz_clear( k_memo->four_k_factorial );
    mpz_clear( k_memo->k_factorial_fourth );
    mpz_clear( k_memo->k_26390 );
    mpz_clear( k_memo->three_ninety_six_to_four_k );
}

void print_memo(struct Memo *k_memo) {
    printf ("Memo { k : %ld || ", k_memo->k );
    gmp_printf ("26390k : %Zd || ", k_memo->k_26390);
    gmp_printf ("396^4k : %Zd || ", k_memo->three_ninety_six_to_four_k);
    gmp_printf ("4k! : %Zd || ", k_memo->four_k_factorial);
    gmp_printf ("k!^4 : %Zd }\n", k_memo->k_factorial_fourth);

}

// (4n)! = (4k!) * 4n * 4n-1 ... 4k+1
void next_four_k_factorial(mpz_t current, unsigned long k, unsigned long n) {
    for (unsigned long i = n; i > k; --i) {
        unsigned long fact_term = 4*i;
        mpz_mul_ui(current, current, fact_term);
        mpz_mul_ui(current, current, fact_term - 1);
        mpz_mul_ui(current, current, fact_term - 2);
        mpz_mul_ui(current, current, fact_term - 3);
    }
}

 // n!^4 = k!^4 * n^4 * (n-1)^4 ... (k+1)^4
void next_k_factorial_fourth(mpz_t current, unsigned long k, unsigned long n) {
    mpz_t n_minus_x_fourth;
    mpz_init( n_minus_x_fourth );
    for (unsigned long i = n; i > k; --i) {
        mpz_ui_pow_ui(n_minus_x_fourth, i, 4UL); // (n-x)^4
        mpz_mul(current, current, n_minus_x_fourth);
    }

    // Free up temp
    mpz_clear( n_minus_x_fourth );
}

// 26390n = 26390k + 26390(n-k)
void next_k_26390(mpz_t current, unsigned long k, unsigned long n) {
    mpz_add_ui(current, current, 26390 * (n-k));
}

// 396^4n = 396^4k * 396^4(n-k)
void next_three_ninety_six_to_four_k(mpz_t current, unsigned long k, unsigned long n) {
    mpz_t three_ninety_six_to_four_n_minus_k;
    mpz_init( three_ninety_six_to_four_n_minus_k );
    mpz_ui_pow_ui(three_ninety_six_to_four_n_minus_k, 396UL, 4UL*(n-k));
    mpz_mul(current, current, three_ninety_six_to_four_n_minus_k);

    // Free up temp
    mpz_clear( three_ninety_six_to_four_n_minus_k );
}

void advance_to_n(struct Memo *k_memo, unsigned long n) {
    next_four_k_factorial(k_memo->four_k_factorial, k_memo->k, n);
    next_k_factorial_fourth(k_memo->k_factorial_fourth, k_memo->k, n);
    next_k_26390(k_memo->k_26390, k_memo->k, n);
    next_three_ninety_six_to_four_k(k_memo->three_ninety_six_to_four_k, k_memo->k, n);
    k_memo->k = n;
}

void update_accumulator(const struct Memo *k_memo, mpfr_t accumulator, int desired_precision) {
    mpfr_t term;
    mpfr_init2( term, desired_precision );
    // Set it to 26390k + 1103
    mpfr_set_z( term, k_memo->k_26390, MPFR_RNDN );
    mpfr_add_ui( term, term, 1103UL, MPFR_RNDN );
    // Divide by 396^4k
    mpfr_div_z( term, term, k_memo->three_ninety_six_to_four_k, MPFR_RNDN );
    // Multiply by 4k!
    mpfr_mul_z( term, term, k_memo->four_k_factorial, MPFR_RNDN );
    // Divide by k!^4
    mpfr_div_z( term, term, k_memo->k_factorial_fourth, MPFR_RNDN );

    // introduce the sqrt(8) /9801 factor
    mpfr_t factor;
    mpfr_init2( factor, desired_precision );
    mpfr_sqrt_ui( factor, 8UL, MPFR_RNDN );
    mpfr_div_ui( factor, factor, 9801UL, MPFR_RNDN );

    // multiply the factor into our term
    mpfr_mul( term, term, factor, MPFR_RNDN );

    // Update accumulator
    mpfr_add( accumulator, accumulator, term, MPFR_RNDN );

    // Free up temp
    mpfr_clear( term );
    mpfr_clear( factor );
}

int check_results(mpfr_t accumulator, long precision) {

    FILE * pi_file = fopen("pi.txt", "r");
    if (pi_file == NULL) {
        printf("can't find pi.txt");
        return EXIT_FAILURE;
    }
    char * buffer;
    char * calc_buffer = NULL;
    long   numbytes;
    // Get the number of bytes
    fseek(pi_file, 0L, SEEK_END);
    numbytes = ftell(pi_file);

    // reset the file position indicator to
    // the beginning of the file
    fseek(pi_file, 0L, SEEK_SET);

    // grab sufficient memory for the
    // buffer to hold the text
    buffer = (char*)calloc(numbytes, sizeof(char));

    // memory error
    if(buffer == NULL)
        return EXIT_FAILURE;

    // printf( "reading file\n" );
    // copy all the text into the buffer
    fread( buffer, sizeof(char), numbytes, pi_file );
    fclose( pi_file );

    // printf( "Getting pi string\n" );
    // Write our calulated pi to the buffer
    long i = 1;
    long * i_ptr = &i;
    calc_buffer = mpfr_get_str( calc_buffer, i_ptr, 10, precision, accumulator, MPFR_RNDN );


    // Check our accuracy
    for (unsigned long overlap = 0; overlap <= precision; ++overlap) {
        // printf("Comparing %c to %c\n", buffer[overlap], calc_buffer[overlap]);
        if ( buffer[overlap] != calc_buffer[overlap] ) {
            printf("Comparing %c to %c\n", buffer[overlap], calc_buffer[overlap]);
            printf("Calculation accurate to %ld decimal places\n", overlap);
            break;
        }
    }

    // free the memory we used for the buffer
    free(buffer);
    mpfr_free_str(calc_buffer);

    return EXIT_SUCCESS;
}

int main( int argc, char *argv[] ) {
    // Handle command line arg for changing precision

    char *char_pointer;
    long precision = DECIMAL_PRECISION;
    if (argc > 1) {
        errno = 0;
        long converted = strtol(argv[1], &char_pointer, 10);

        if (errno == 0 && *char_pointer == '\0' && converted > 0 && converted <= __LONG_MAX__) {
            precision = converted;
        }
    }

    // How many bits are needed to achieve the desired decimal precision?
    const int desired_precision = (int)ceil( precision * ( log( 10 ) / log( 2 ) ) );
    mpfr_t accumulator;
    struct Memo k_memo;

    mpfr_init2( accumulator, desired_precision );
    mpfr_set_d( accumulator, 0.0, MPFR_RNDN );
    memo_init( &k_memo );

    unsigned long iterations = precision / 7;
    printf("Making %ld iterations...\n", iterations);
    // print_memo(&k_memo);
    clock_t start_time = clock();
    clock_t checkpoint;
    for (unsigned long k = 0; k <= iterations; ++k) {
        if (k == iterations || (k > 0 && k % 10000 == 0)) {
            checkpoint = clock() - start_time;
            printf("%ld iterations complete in %fs\n", k, (double)checkpoint / CLOCKS_PER_SEC);
        }
        advance_to_n( &k_memo, k );
        // print_memo(&k_memo);
        update_accumulator( &k_memo, accumulator, desired_precision);
    }

    // account for the 1/pi thing
    mpfr_ui_div( accumulator, 1UL, accumulator, MPFR_RNDN );

    // Print the answer.
    printf( "pi = " );
    mpfr_out_str( stdout, 10, 50, accumulator, MPFR_RNDN );
    printf( "\n" );

    int result = check_results(accumulator, precision);

    // Release resources associated with the multi-precision floats.
    mpfr_clear( accumulator );
    // Clear the memo struct
    memo_destroy(&k_memo);

    return result;
}
