/*
 * \file sum.c
 * \brief An example of various ways to sum an array of double.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <locale.h>

#include <pthread.h>
#if defined(__GLIBC__) || defined(__CYGWIN__)
#include <sys/sysinfo.h>
#endif

#include "Timer.h"

#define SIZE            100000000    // The size of the array to sum.
#define ITERATIONS      100          // The number of times we'll do the summation for timing purposes.
#define TOTAL_VOLUME    1000000000   // How much to churn

// -----------
// Serial Sums
// -----------

double sum_simple( const double *array, size_t size )
{
    double accumulator = 0.0;
    for( size_t i = 0; i < size; ++i ) {
        accumulator += array[i];
    }
    return accumulator;
}

double sum_recursive( const double *array, size_t size )
{
    int midpoint = size / 2;

    if( size == 1 ) return array[0];
    return sum_recursive( array, midpoint ) +
           sum_recursive( array + midpoint, size - midpoint );
}

double sum_hybrid( const double *array, size_t size )
{
    size_t midpoint = size / 2;

    if( size <= 10000 ) return sum_simple( array, size );
    return sum_hybrid( array, midpoint ) +
           sum_hybrid( array + midpoint, size - midpoint );
}

// ------------
// Parallel Sum
// ------------

// Structure to define the range of the array processed by a single thread.
struct PThreadWorkUnit {
    const double *start;
    size_t size;
};


// The thread function that adds up all elements in a given range of a given array.
void *summer( void *arg )
{
    struct PThreadWorkUnit *range = (struct PThreadWorkUnit *)arg;
    double accumulator = 0.0;
    double *result;

    for( size_t i = 0; i < range->size; ++i ) {
        accumulator += range->start[i];
    }

    result = (double *)malloc( sizeof(double) );
    *result = accumulator;
    return result;
}


//! Add all the integers in the array in parallel.
/*!
 * This function uses exactly two threads to add the elements in a (large) array. There is no
 * attempt here to scale the number of threads to match the underlying hardware.
 *
 * \param array A pointer to the array (or section of an array) to sum.
 * \param size The number of elements in the array to sum.
 * \return The sum of all the array elements. There is no checking for overflow.
 */
double sum_parallel( const double *array, size_t size )
{
    struct PThreadWorkUnit ranges[2];   // Define the subproblems (thread arguments)
    pthread_t threads[2];          // Thread handles.
    double    accumulator = 0.0;   // The final answer.
    void     *result;              // Used to hold the values returned by the threads.

    // Split the problem.
    ranges[0].start = array;
    ranges[0].size  = size / 2;
    ranges[1].start = array + size / 2;
    ranges[1].size  = size - size / 2;

    // Start the worker threads.
    pthread_create( &threads[0], NULL, summer, &ranges[0] );
    pthread_create( &threads[1], NULL, summer, &ranges[1] );

    // Gather the results.
    pthread_join( threads[0], &result );
    accumulator += *(double*)result;
    free( result );

    pthread_join( threads[1], &result );
    accumulator += *(double*)result;
    free( result );

    // Return aggregated result.
    return accumulator;
}


//! Add all the integers in the array in parallel using a dynamic number of processors.
/*!
 * This version scales the number of threads to account for the capabilities of the underlying
 * hardware.
 *
 * \param array A pointer to the array (or section of an array) to sum.
 * \param size The number of elements in the array to sum.
 * \return The sum of all the array elements. There is no checking for overflow.
 */
double sum_dynamic( const double *array, size_t size, size_t processor_count )
{

    struct PThreadWorkUnit *ranges =
        (struct PThreadWorkUnit *)malloc( processor_count * sizeof(struct PThreadWorkUnit) );
    pthread_t *threads =
        (pthread_t *)malloc( processor_count * sizeof(pthread_t) );

    size_t  i;
    size_t  chunk_size;
    double  accumulator = 0;
    void   *result;

    // Split the problem.
    chunk_size = size / processor_count;
    for( i = 0; i < processor_count; ++i ) {
        ranges[i].start = array + ( i * chunk_size );
        ranges[i].size = chunk_size;
    }
    // The following line assigns the remainder elements to the last thread.
    ranges[processor_count - 1].size = size - (processor_count - 1) * chunk_size;

    // Start the worker threads.
    for( i = 0; i < processor_count; ++i ) {
        pthread_create( &threads[i], NULL, summer, &ranges[i] );
    }

    // Gather the results.
    for( i = 0; i < processor_count; ++i ) {
        pthread_join( threads[i], &result );
        accumulator += *(double *)result;
        free( result );
    }

    // Release dynamic memory.
    free( threads );
    free( ranges );

    // Return aggregated result.
    return accumulator;
}

// -------
// Testing
// -------

// The runner function exercises either the serial or the parallel version.
double runner( const char *tag, const double *array, size_t size, size_t iterations, double ( *function )( const double *, size_t, size_t ), size_t processor_count )
{
    double sum;
    Timer  stopwatch;
    double seconds;

    Timer_initialize( &stopwatch );
    Timer_start( &stopwatch );
    for( int i = 0; i < iterations; ++i ) {
        sum = function( array, size, processor_count );
    }
    Timer_stop( &stopwatch );
    seconds = (double)Timer_time( &stopwatch ) / 1000.0;
    seconds /= iterations;
    // printf( "Sum (%s) = %f (%f seconds)\n", tag, sum, seconds );

    return seconds;
}


int main( int argc, char *argv[] )
{
    // Handle command line arg for changing size

    char *char_pointer;
    long size = SIZE;
    if (argc > 1) {
        errno = 0;
        long converted = strtol(argv[1], &char_pointer, 10);

        if (errno == 0 && *char_pointer == '\0' && converted > 0 && converted <= __LONG_MAX__) {
            size = converted;
        }
    }
    #if defined(__GLIBC__) || defined(__CYGWIN__)
    // A glibc-specific function.
    int processor_count = get_nprocs( );
    #else
    // A POSIX threads extension (not supported by glibc, but available, e.g., on Windows).
    int processor_count = pthread_num_processors_np( );
    #endif

    long iterations = TOTAL_VOLUME / size;

    int rc = EXIT_SUCCESS;
    double *p = (double *)malloc( size * sizeof(double) );

    /*
    #if defined(__GLIBC__) || defined(__CYGWIN__)
    // A glibc-specific function.
    int processor_count = get_nprocs( );
    #else
    // A POSIX threads extension (not supported by glibc, but available, e.g., on Windows).
    int processor_count = pthread_num_processors_np( );
    #endif

    printf( "This environment has %d processing elements.\n", processor_count );
    */

    if( p == NULL ) {
        printf( "Unable to allocate memory for array!\n" );
        rc = EXIT_FAILURE;
    }
    else {
        Timer  stopwatch;
        double seconds;

        Timer_initialize( &stopwatch );
        Timer_start( &stopwatch );
        // Fill array with known data.
        for( size_t i = 0; i < size; ++i ) {
            p[i] = 1.0;
        }
        Timer_stop( &stopwatch );
        seconds = (double)Timer_time( &stopwatch ) / 1000.0;
        printf( "Filled Array Size %ld in %lf seconds\n", size, seconds );

        // Try the various ways of summing it.
        // runner( "simp", p, size, iterations, sum_simple );
        // runner( "recu", p, size, iterations, sum_recursive );
        // double hybrid_time = runner( "hybr", p, size, iterations, sum_hybrid );
        // double parallel_time = runner( "para", p, size, iterations, sum_parallel );

	double single_threaded_time = 0;

        for (size_t i = 1; i <= processor_count; ++i) {
            double dynamic_time = runner( "dyna", p, size, iterations, sum_dynamic, i );
            double dynamic_time_two = runner( "dyna", p, size, iterations, sum_dynamic, i );
            double dynamic_time_three = runner( "dyna", p, size, iterations, sum_dynamic, i );

            // Grab the best of 3 runs
            if (dynamic_time_two <= dynamic_time && dynamic_time_two <= dynamic_time_three) {
                dynamic_time = dynamic_time_two;
            } else if (dynamic_time_three <= dynamic_time && dynamic_time_three <= dynamic_time_two) {
                dynamic_time = dynamic_time_three;
            }

            if (i == 1) {
                single_threaded_time = dynamic_time;
	            printf("Array Size: %ld  Threads: %ld  Time: %lf\n", size, i, dynamic_time);
            } else {
	            printf("Array Size: %ld  Threads: %ld  Time: %lf  Speedup: %lf\n", size, i, dynamic_time, single_threaded_time / dynamic_time);
            }
	}
        // printf( "%16ld,%16lf,%16lf,%16lf,%16lf,%16lf\n", size, hybrid_time, parallel_time, dynamic_time, hybrid_time / parallel_time, hybrid_time / dynamic_time);
        free( p );
    }

    return rc;
}
