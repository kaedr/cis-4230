/*!
 * \file   gaussian.h
 * \brief  Interface to a Gaussian Elimination solver.
 * \author (C) Copyright 2024 by Peter Chapin <pchapin@vermontstate.edu>
 */

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <stdlib.h>

// The symbol __STDC_NO_VLA__ is part of the C 2011 standard.
#ifdef __STDC_NO_VLA___
#error C99-style variable length arrays are required, but are not supported by this compiler.
#endif

// Change this type alias to change the data type of the matrix elements.
typedef double floating_type;

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
enum GaussianResult gaussian_solve( size_t size, floating_type (* restrict a)[size], floating_type * restrict b, int selection );

#endif
