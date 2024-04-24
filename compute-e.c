/*!
 * \file    main.c
 * \brief   Sample program that illustrates how to use GMP and MPFR.
 * \author  Peter Chapin <peter.chapin@vermontstate.edu>
 *
 * This program computes the value of `e`, the base of natural logarithms. It uses the following
 * infinite series expansion of e^x:
 *
 * e^x = \sum_{k = 0}^{\infty} \frac{x^k}{k!}
 *
 * Setting `x = 1` gives:
 *
 * e = \sum_{k = 0}^{\infty} \frac{1}{k!} = 1 + 1 + 1/2 + 1/6 + 1/24 + 1/120 + ... = 2.71828...
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

// The maximum number of terms in the series to add.
#define MAX_TERM_INDEX    500
#define DECIMAL_PRECISION 100

int main( void )
{
    int k;
    mpz_t  k_factorial;
    mpfr_t accumulator;
    mpfr_t term;

    // The `mpz_init` function sets the integer to zero. We initialize it to 0!.
    // The function `mpz_init_set_ui` can be used to combine these operations.
    mpz_init( k_factorial );
    mpz_set_ui( k_factorial, 1UL );

    // How many bits are needed to achieve the desired decimal precision?
    const int desired_precision = (int)ceil( DECIMAL_PRECISION * ( log( 10 ) / log( 2 ) ) );

    // The mpfr_init2 function sets the floating point value to NaN. We initialize it to 1.0.
    // As with GMP, there are also functions that initialize the number and set it to a value.
    mpfr_init2( accumulator, desired_precision );
    mpfr_set_d( accumulator, 1.0, MPFR_RNDN );
    mpfr_init2( term, desired_precision );

    // Compute the result. The k = 0 term is already handled (accumulator initialized to 1.0).
    for( k = 1; k <= MAX_TERM_INDEX; ++k ) {
        mpz_mul_ui( k_factorial, k_factorial, k );             // Update k!.
        mpfr_set_d( term, 1.0, MPFR_RNDN );                    // The numerator.
        mpfr_div_z( term, term, k_factorial, MPFR_RNDN );      // Compute 1/k!.
        mpfr_add( accumulator, accumulator, term, MPFR_RNDN ); // Add term to accumulator.
    }

    // Print the answer.
    printf( "e = " );
    mpfr_out_str( stdout, 10, 0, accumulator, MPFR_RNDN );
    printf( "\n" );

    // Release resources associated with the multi-precision floats.
    mpfr_clear( term );
    mpfr_clear( accumulator );

    // Release resources associated with the multi-precision integers.
    mpz_clear( k_factorial );
    return EXIT_SUCCESS;
}
