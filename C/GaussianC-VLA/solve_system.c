/*!
 *  \file   solve_system.c
 *  \brief  Solve a large system of simultaneous equations.
 *  \author (C) Copyright 2024 by Peter Chapin <pchapin@vermontstate.edu>
 */

#include <stdlib.h>
#include <stdio.h>

#include "gaussian.h"
#include "Timer.h"


int main( int argc, char *argv[] )
{
    FILE   *input_file;
    size_t  size;

    if( argc != 2 ) {
        printf( "Error: Expected the name of a system definition file.\n" );
        return EXIT_FAILURE;
    }

    // Open the file.
    if( (input_file = fopen( argv[1], "r" )) == NULL ) {
        printf("Error: Can not open the system definition file.\n");
        return EXIT_FAILURE;
    }

    // Get the size.
    fscanf( input_file, "%zu", &size );

    // Allocate the arrays on the stack... except this overflows the stack for large systems.
    //floating_type a[size][size];
    //floating_type b[size];

    // Allocate the arrays dynamically.
    typedef floating_type row_t[size];
    row_t *a = (row_t *)malloc( size * size * sizeof( floating_type ) );
    floating_type *b = (floating_type *)malloc( size * sizeof( floating_type ) );

    // Get coefficients.
    // Note that the format specifier used here, `%lf`, assumes the matrix elements have type
    // double. See the declaration of `floating_type` at the top of gaussian.h.
    //
    for( size_t i = 0; i < size; ++i ) {
        for( size_t j = 0; j < size; ++j ) {
            fscanf( input_file, "%lf", &a[i][j] );
        }
        fscanf( input_file, "%lf", &b[i] );
    }
    fclose( input_file );
    //printf( "\nFinished reading %s\n", argv[1] );

    // Do the calculations.
    Timer stopwatch;
    Timer_initialize( &stopwatch );
    Timer_start( &stopwatch );
    enum GaussianResult result = gaussian_solve( size, a, b );
    Timer_stop( &stopwatch );

    // Display the results.
    switch( result ) {
    case gaussian_success:
        // printf( "\nSolution is\n" );
        // for( size_t i = 0; i < size; ++i ) {
        //     printf( " x[%4zu] = %9.5f\n", i, b[i] );
        // }
        printf( "\nExecution time = %ld milliseconds\n", Timer_time( &stopwatch ) );
        break;

    case gaussian_error:
        printf( "Parameter problem in call to gaussian_solve( )\n" );
        break;

    case gaussian_degenerate:
        printf( "System is degenerate. It does not have a unique solution.\n" );
        break;
    }

    // Clean up the dynamically allocated space.
    free( a );
    free( b );
    return EXIT_SUCCESS;
}
