/*!
 * \file   gaussian.c
 * \brief  A Gaussian Elimination solver.
 * \author (C) Copyright 2024 by Peter Chapin <pchapin@vermontstate.edu>
 *
 * This is the serial version of the algorithm.
 */

#include <math.h>
#include <string.h>

#include "gaussian.h"

// For profiling, it is best for all functions to be public.
#define PRIVATE // static
#define PUBLIC

//! Does the elimination step of reducing the system. O(n^3)
PRIVATE enum GaussianResult elimination( size_t size, floating_type (* restrict a)[size], floating_type * restrict b )
{
    //floating_type *temp_array = (floating_type *)malloc( size * sizeof(floating_type) );
    floating_type  temp_array[size];
    size_t         i, j, k;
    floating_type  temp, m;

    for( i = 0; i < size - 1; ++i ) {

        // Find the row with the largest value of |a[j][i]|, j = i, ..., n - 1
        k = i;
        m = fabs( a[i][i] );
        for( j = i + 1; j < size; ++j ) {
            if( fabs( a[j][i] ) > m ) {
                k = j;
                m = fabs( a[j][i] );
            }
        }

        // Check for |a[k][i]| zero.
        // TODO: The value 1.0E-6 is arbitrary. A more disciplined value should be used.
        if( fabs( a[k][i] ) <= 1.0E-6 ) {
            //free( temp_array );
            return gaussian_degenerate;
        }

        // Exchange row i and row k, if necessary.
        if( k != i ) {
            memcpy( temp_array, a[i], size * sizeof( floating_type ) );
            memcpy( a[i], a[k], size * sizeof( floating_type ) );
            memcpy( a[k], temp_array, size * sizeof( floating_type ) );
            
            // Exchange corresponding elements of b.
            temp = b[i];
            b[i] = b[k];
            b[k] = temp;
        }

        // Subtract multiples of row i from subsequent rows.
        for( j = i + 1; j < size; ++j ) {
            m = a[j][i] / a[i][i];
            for( k = 0; k < size; ++k )
                a[j][k] -= m * a[i][k];
            b[j] -= m * b[i];
        }
    }
    //free( temp_array );
    return gaussian_success;
}


//! Does the back substitution step of solving the system. O(n^2)
PRIVATE enum GaussianResult back_substitution( size_t size, floating_type (* restrict a)[size], floating_type * restrict b )
{
    floating_type sum;
    size_t        i, j;
    size_t        counter;

    // We can't count i down from size - 1 to zero (inclusive) because it is unsigned.
    for( counter = 0; counter < size; ++counter ) {
        i = ( size - 1 ) - counter;
        // TODO: The value 1.0E-6 is arbitrary. A more disciplined value should be used.
        if( fabs( a[i][i] ) <= 1.0E-6 ) {
            return gaussian_degenerate;
        }

        sum = b[i];
        for( j = i + 1; j < size; ++j ) {
            sum -= a[i][j] * b[j];
        }
        b[i] = sum / a[i][i];
    }
    return gaussian_success;
}


PUBLIC enum GaussianResult gaussian_solve( size_t size, floating_type (* restrict a)[size], floating_type * restrict b )
{
    // We can deal with a 1x1 system, but not an empty system.
    if( size == 0 ) return gaussian_error;

    enum GaussianResult return_code = elimination( size, a, b );
    if( return_code == gaussian_success )
        return_code = back_substitution( size, a, b );
    return return_code;
}
