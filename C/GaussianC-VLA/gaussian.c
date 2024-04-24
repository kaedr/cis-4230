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

#include "ThreadPool.h"
#include "gaussian.h"

// For profiling, it is best for all functions to be public.
#define PRIVATE // static
#define PUBLIC

#define PROCESSOR_COUNT 8

//! Does the elimination step of reducing the system. O(n^3)
PRIVATE enum GaussianResult serial_elimination( size_t size, floating_type (* restrict a)[size], floating_type * restrict b )
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
        if (i % 2 == 0) {
            for( j = i + 1; j < size; ++j ) {
                m = a[j][i] / a[i][i];
                for( k = 0; k < size; ++k ) {
                    a[j][k] -= m * a[i][k];
                }
                b[j] -= m * b[i];
            }
        } else {
            for( j = size - 1; j > i; --j ) {
                m = a[j][i] / a[i][i];
                for( k = 0; k < size; ++k ) {
                    a[j][k] -= m * a[i][k];
                }
                b[j] -= m * b[i];
            }
        }
    }
    //free( temp_array );
    return gaussian_success;
}

// Structure to define the data processed by a single thread.
struct PThreadWorkUnit {
    floating_type *a;
    floating_type *b;
    size_t current;
    size_t start;
    size_t stop;
    size_t size;
};

void * p_thread_chunk_elimination( void *arg ) {
    struct PThreadWorkUnit *unit = (struct PThreadWorkUnit *)arg;

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

enum GaussianResult p_thread_elimination( size_t size, floating_type (* restrict a)[size], floating_type * restrict b ) {

    int processor_count = PROCESSOR_COUNT;
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

        struct PThreadWorkUnit *ranges =
            (struct PThreadWorkUnit *)malloc( processor_count * sizeof(struct PThreadWorkUnit) );
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
            pthread_create( &threads[h], NULL, p_thread_chunk_elimination, &ranges[h] );
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


// Structure to define the data processed by a single thread.
struct BarrierWorkUnit {
    floating_type *a;
    floating_type *b;
    floating_type *temp_array;
    size_t size;
    size_t offset;
    enum GaussianResult result;
};

pthread_barrier_t iteration_barrier;
pthread_barrier_t work_barrier;

void * barrier_work( void *arg) {
    struct BarrierWorkUnit *unit = (struct BarrierWorkUnit *)arg;

    const size_t size = unit->size;
    const size_t offset = unit->offset;

    floating_type (* restrict a)[size] = unit->a;
    floating_type * restrict b = unit->b;
    floating_type (* restrict temp_array)[size] = unit->temp_array;

    size_t         start, stop;
    size_t         problem_size, chunk_size;
    size_t         i, j, k;
    floating_type  temp, m;

    for( i = 0; i < size - 1; ++i ) {
        if ( pthread_barrier_wait( &iteration_barrier ) == PTHREAD_BARRIER_SERIAL_THREAD ) {
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
                unit-> result = gaussian_degenerate;
                return NULL;
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

        }

        pthread_barrier_wait( &work_barrier );

        problem_size = ( size - ( i + 1 ));
        chunk_size = problem_size / PROCESSOR_COUNT;
        start = i + 1 + offset * chunk_size;
        if (offset == PROCESSOR_COUNT - 1) {
            stop = size;
        } else {
            stop = start + chunk_size;
        }
        // Subtract multiples of row i from subsequent rows.
        if (i % 2 == 0) {
            for( j = start; j < stop; ++j ) {
                m = a[j][i] / a[i][i];
                for( k = 0; k < size; ++k ) {
                    a[j][k] -= m * a[i][k];
                }
                b[j] -= m * b[i];
            }
        } else {
            for( j = stop - 1; j >= start; --j ) {
                m = a[j][i] / a[i][i];
                for( k = 0; k < size; ++k ) {
                    a[j][k] -= m * a[i][k];
                }
                b[j] -= m * b[i];
            }
        }
    }

    unit-> result = gaussian_success;
    return NULL;
}

//! Does the elimination step of reducing the system. O(n^3)
PRIVATE enum GaussianResult barrier_elimination( size_t size, floating_type (* restrict a)[size], floating_type * restrict b )
{
    //floating_type *temp_array = (floating_type *)malloc( size * sizeof(floating_type) );
    floating_type  temp_array[size];

    struct BarrierWorkUnit *units =
            (struct BarrierWorkUnit *)malloc( PROCESSOR_COUNT * sizeof(struct BarrierWorkUnit) );
    pthread_t *threads =
        (pthread_t *)malloc( PROCESSOR_COUNT * sizeof(pthread_t) );

    pthread_barrier_init( &iteration_barrier, NULL, PROCESSOR_COUNT );
    pthread_barrier_init( &work_barrier, NULL, PROCESSOR_COUNT );

    // Create a thread for each CPU and set it working on its work unit.
    for( int offset = 0; offset < PROCESSOR_COUNT; ++offset ) {
        units[offset].a = a;
        units[offset].b = b;
        units[offset].temp_array = temp_array;
        units[offset].size = size;
        units[offset].offset = offset;

        pthread_create( &threads[offset], NULL, barrier_work, &units[offset]);
    }


    for( int h = 0; h < PROCESSOR_COUNT; ++h ) {
        pthread_join( threads[h], NULL );
        if (units[h].result == gaussian_degenerate) {
            return gaussian_degenerate;
        }
    }

    // Release dynamic memory.
    free( threads );
    free( units );
    return gaussian_success;
}


// Structure to define the data processed by a single thread.
struct PoolWorkUnit {
    floating_type *a;
    floating_type *b;
    size_t current;
    size_t start;
    size_t stop;
    size_t size;
};

void * pool_chunk_elimination( void *arg ) {
    struct PoolWorkUnit *unit = (struct PoolWorkUnit *)arg;

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

enum GaussianResult pool_elimination( size_t size, floating_type (* restrict a)[size], floating_type * restrict b ) {

    int processor_count = 8;
    floating_type  temp_array[size];
    size_t         i, j, k;
    floating_type  temp, m;
    size_t  chunk_size;

    ThreadPool pool;
    ThreadPool_initialize(&pool);



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

        struct PoolWorkUnit *ranges =
            (struct PoolWorkUnit *)malloc( processor_count * sizeof(struct PoolWorkUnit) );
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
            threads[h] = ThreadPool_start( &pool, pool_chunk_elimination, &ranges[h]);
        }

        for( int h = 0; h < processor_count; ++h ) {
            ThreadPool_result(&pool, threads[h]);
        }

        // Release dynamic memory.
        free( threads );
        free( ranges );

    }

    ThreadPool_destroy(&pool);

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


PUBLIC enum GaussianResult gaussian_solve( size_t size, floating_type (* restrict a)[size], floating_type * restrict b, int selection )
{
    // We can deal with a 1x1 system, but not an empty system.
    if( size == 0 ) return gaussian_error;
    enum GaussianResult return_code;

    switch (selection)
    {
    // Serial
    case 1:
    case '1':
        return_code = serial_elimination( size, a, b );
        break;
    // p_thread
    case 2:
    case '2':
        return_code = p_thread_elimination( size, a, b );
        break;
    // Barrier
    case 3:
    case '3':
        return_code = barrier_elimination( size, a, b );
        break;
    // Thread Pool
    case 4:
    case '4':
        return_code = pool_elimination( size, a, b );
        break;

    default:
        return gaussian_error;
    }

    if( return_code == gaussian_success )
        return_code = back_substitution( size, a, b );
    return return_code;
}
