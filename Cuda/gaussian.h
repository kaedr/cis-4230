/*!
 * \file   gaussian.h
 * \brief  Interface to a Gaussian Elimination solver.
 * \author (C) Copyright 2024 by Peter Chapin <pchapin@vermontstate.edu>
 */

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <stdlib.h>

// Change this type alias to change the data type of the matrix elements.
typedef double floating_type;

// Macros for handling matrixes.
//
// These macros manipulate a linear array as if it was a two dimensional array. Note that C99
// variable-length arrays can be used to deal with this problem much more elegantly.
//
// TODO: Create a matrix abstraction? Or would the overhead of doing so be too great?
#define MATRIX_MAKE( size )  ((floating_type *)malloc( (size) * (size) * sizeof( floating_type ) ))
#define MATRIX_DESTROY( matrix )                       ( free( matrix ) )
#define MATRIX_GET( matrix, size, row, column )        ( (matrix)[(row)*(size) + (column)] )
#define MATRIX_GET_REF( matrix, size, row, column )    (&(matrix)[(row)*(size) + (column)] )
#define MATRIX_GET_ROW( matrix, size, row )            (&(matrix)[(row)*(size)] )
#define MATRIX_PUT( matrix, size, row, column, value ) ( (matrix)[(row)*(size) + (column)] = (value) )

enum GaussianResult {
    gaussian_success,     // The system was solved normally.
    gaussian_error,       // A problem with the parameters was detected.
    gaussian_degenerate   // The system is degenerate and does not have a unique solution.
};

//! Gaussian Elimination solver.
/*!
 * \param a A pointer to the matrix of coefficients in row-major order.
 * \param b A pointer to the driving vector.
 * \returns gaussian_success if the system is solved.
 *
 * This function solves the system in place. If it is successful, the driving vector is replaced
 * with the solution. If this function is not successful, the matrix of coefficients and the
 * driving vector may be in a partially modified state.
 */
enum GaussianResult gaussian_solve( size_t size, floating_type *a, floating_type *b );

#endif
