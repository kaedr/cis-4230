/*!
 * \file   gaussian.c
 * \brief  A Gaussian Elimination solver.
 * \author (C) Copyright 2024 by Peter Chapin <pchapin@vermontstate.edu>
 *
 * This is the serial version of the algorithm.
 */

#include <math.h>
#include <string.h>
#include <pthread.h>
#include <stdio.h>

#include "gaussian.h"

// For profiling, it is best for all functions to be public.
#define PRIVATE // static
#define PUBLIC

// Structure to define the data processed by a single thread.
struct WorkUnit {
    floating_type *a;
    floating_type *b;
    size_t current;
    size_t start;
    size_t stop;
    size_t size;
};

void * chunk_elimination( void *arg ) {
    struct WorkUnit *unit = (struct WorkUnit *)arg;

    const size_t size = unit->size;
    const size_t current = unit->current;
    const size_t start = unit->start;
    const size_t stop = unit->stop;

    floating_type (* restrict a)[size] = unit->a;
    floating_type * restrict b = unit->b;
    size_t         j, k;
    floating_type  m;

    // Subtract multiples of row i from subsequent rows.
    if (start % 2 == 0) {
        for( j = start; j < stop; ++j ) {
            m = a[j][current] / a[current][current];
            for( k = 0; k < size; ++k ) {
                a[j][k] -= m * a[current][k];
            }
            b[j] -= m * b[current];
        }
    } else {
        for( j = stop - 1; j >= start; --j ) {
            m = a[j][current] / a[current][current];
            for( k = 0; k < size; ++k ) {
                a[j][k] -= m * a[current][k];
            }
            b[j] -= m * b[current];
        }

    }

    return NULL;
}

enum GaussianResult parallel_elimination( size_t size, floating_type (* restrict a)[size], floating_type * restrict b ) {

    int processor_count = 8;
    floating_type  temp_array[size];
    size_t         i, j, k;
    floating_type  temp, m;
    size_t  chunk_size;

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

        struct WorkUnit *ranges =
            (struct WorkUnit *)malloc( processor_count * sizeof(struct WorkUnit) );
        pthread_t *threads =
            (pthread_t *)malloc( processor_count * sizeof(pthread_t) );

        // Split the problem.
        size_t problem_size = ( size - ( i + 1 ));
        chunk_size = problem_size / processor_count;
        for( size_t x = 0; x < processor_count; ++x ) {
            ranges[x].a = a;
            ranges[x].b = b;
            ranges[x].start = i + 1 + x * chunk_size;
            ranges[x].stop = ranges[x].start + chunk_size;
            ranges[x].current = i;
            ranges[x].size = size;
        }
        // The following line assigns the remainder elements to the last thread.
        ranges[processor_count - 1].stop = size;

        // Start the worker threads.
        for( int h = 0; h < processor_count; ++h ) {
            pthread_create( &threads[h], NULL, chunk_elimination, &ranges[h] );
        }

        for( int h = 0; h < processor_count; ++h ) {
            pthread_join( threads[h], NULL );
        }

        // Release dynamic memory.
        free( threads );
        free( ranges );

    }

    return gaussian_success;
}

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
            for( k = 0; k < size; ++k ) {
                a[j][k] -= m * a[i][k];
            }
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

    enum GaussianResult return_code = parallel_elimination( size, a, b );
    if( return_code == gaussian_success )
        return_code = back_substitution( size, a, b );
    return return_code;
}
