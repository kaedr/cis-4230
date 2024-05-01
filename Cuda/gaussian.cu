/*!
 * \file   gaussian.c
 * \brief  A Gaussian Elimination solver.
 * \author (C) Copyright 2024 by Peter Chapin <pchapin@vermontstate.edu>
 *
 * This is the serial version of the algorithm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda_runtime.h>

#include "gaussian.h"

#define PRIVATE static
#define PUBLIC

__global__ void elimination_kernel( size_t size, floating_type *a, floating_type *b, size_t i ) {
    size_t         k;
    floating_type  m;

    int my_id = threadIdx.x;
    int j = i + 1 + my_id;
    // if (j >= size) {
    //     printf("Ignoring j=%d\n", j);
    //     return;
    // }
    // Synchronize
    // __syncthreads();
    // Subtract multiples of row i from subsequent rows.
    // for( j = i + 1; j < size; ++j ) {
    //     if (!j % my_id == 0) {
    //         continue;
    //     }
    //     m = MATRIX_GET( a, size, j, i ) / MATRIX_GET( a, size, i, i );
    //     for( k = 0; k < size; ++k )
    //         MATRIX_PUT( a, size, j, k, MATRIX_GET( a, size, j, k ) - m * MATRIX_GET( a, size, i, k ) );
    //     b[j] -= m * b[i];
    // }
    // printf("i: %ld, j: %d my_id: %d\n", i, j, my_id);
    m = MATRIX_GET( a, size, j, i ) / MATRIX_GET( a, size, i, i );
    for( k = 0; k < size; ++k )
        MATRIX_PUT( a, size, j, k, MATRIX_GET( a, size, j, k ) - m * MATRIX_GET( a, size, i, k ) );
    b[j] -= m * b[i];

}

//! Does the elimination step of reducing the system. O(n^3)
PRIVATE enum GaussianResult elimination( size_t size, floating_type *a, floating_type *b )
{
    floating_type *temp_array = (floating_type *)malloc( size * sizeof(floating_type) );
    size_t         i, j, k;
    floating_type  temp, m;

    floating_type *dev_a = MATRIX_MAKE( size );
    floating_type *dev_b = (floating_type *)malloc( size * sizeof( floating_type ) );
    cudaMalloc( (void **)&dev_a, size * size * sizeof(double) );
    cudaMalloc( (void **)&dev_b, size * sizeof(double) );

    for( i = 0; i < size - 1; ++i ) {

        // Find the row with the largest value of |a[j][i]|, j = i, ..., n - 1
        k = i;
        m = fabs( MATRIX_GET( a, size, i, i ) );
        for( j = i + 1; j < size; ++j ) {
            if( fabs( MATRIX_GET( a, size, j, i ) ) > m ) {
                k = j;
                m = fabs( MATRIX_GET( a, size, j, i ) );
            }
        }

        // Check for |a[k][i]| zero.
        // TODO: The value 1.0E-6 is arbitrary. A more disciplined value should be used.
        if( fabs( MATRIX_GET( a, size, k, i ) ) <= 1.0E-6 ) {
            free( temp_array );
            return gaussian_degenerate;
        }

        // Exchange row i and row k, if necessary.
        if( k != i ) {
            memcpy( temp_array, MATRIX_GET_ROW( a, size, i ), size * sizeof( floating_type ) );
            memcpy( MATRIX_GET_ROW( a, size, i ), MATRIX_GET_ROW( a, size, k ), size * sizeof( floating_type ) );
            memcpy( MATRIX_GET_ROW( a, size, k ), temp_array, size * sizeof( floating_type ) );

            // Exchange corresponding elements of b.
            temp = b[i];
            b[i] = b[k];
            b[k] = temp;
        }

        // Subtract multiples of row i from subsequent rows.
        // for( j = i + 1; j < size; ++j ) {
        //     m = MATRIX_GET( a, size, j, i ) / MATRIX_GET( a, size, i, i );
        //     for( k = 0; k < size; ++k )
        //         MATRIX_PUT( a, size, j, k, MATRIX_GET( a, size, j, k ) - m * MATRIX_GET( a, size, i, k ) );
        //     b[j] -= m * b[i];
        // }

        cudaMemcpy( dev_a, a, size * size * sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy( dev_b, b, size * sizeof(double), cudaMemcpyHostToDevice );
        elimination_kernel<<<1, size - 1 - i>>>( size, dev_a, dev_b, i );
        cudaMemcpy( a, dev_a, size * size * sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy( b, dev_b, size * sizeof(double), cudaMemcpyDeviceToHost );
    }
    free( temp_array );
    cudaFree( dev_a );
    cudaFree( dev_b );
    return gaussian_success;
}

//! Does the back substitution step of solving the system. O(n^2)
PRIVATE enum GaussianResult back_substitution( size_t size, floating_type *a, floating_type *b )
{
    floating_type sum;
    size_t        i, j;
    size_t        counter;

    // We can't count i down from size - 1 to zero (inclusive) because it is unsigned.
    for( counter = 0; counter < size; ++counter ) {
        i = ( size - 1 ) - counter;
        // TODO: The value 1.0E-6 is arbitrary. A more disciplined value should be used.
        if( fabs( MATRIX_GET( a, size, i, i ) ) <= 1.0E-6 ) {
            return gaussian_degenerate;
        }

        sum = b[i];
        for( j = i + 1; j < size; ++j ) {
            sum -= MATRIX_GET( a, size, i, j ) * b[j];
        }
        b[i] = sum / MATRIX_GET( a, size, i, i );
    }
    return gaussian_success;
}


PUBLIC enum GaussianResult gaussian_solve( size_t size, floating_type *a, floating_type *b )
{
    // We can deal with a 1x1 system, but not an empty system.
    if( size == 0 ) return gaussian_error;

    enum GaussianResult return_code = elimination( size, a, b );
    if( return_code == gaussian_success )
        return_code = back_substitution( size, a, b );
    return return_code;
}
