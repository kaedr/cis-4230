
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>

#include "Timer.h"



// Change this type alias to change the data type of the matrix elements.
typedef double floating_type;

enum GaussianResult {
    gaussian_success,     // The system was solved normally.
    gaussian_error,       // A problem with the parameters was detected.
    gaussian_degenerate   // The system is degenerate and does not have a unique solution.
};

enum GaussianResult elimination( size_t size, floating_type (* restrict a)[size], floating_type * restrict b )
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

//! Does the back substitution step of solving the system. O(n^2)
enum GaussianResult back_substitution( size_t size, floating_type (* __restrict__ a)[size], floating_type * __restrict__ b )
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

enum GaussianResult gaussian_solve( size_t size, floating_type (* __restrict__ a)[size], floating_type * __restrict__ b )
{
    // We can deal with a 1x1 system, but not an empty system.
    if( size == 0 ) return gaussian_error;

    cl_uint platform_count = 0;

    // Compute the size of the data.
    size_t datasize = sizeof(int) * size;

    // Use this to check the output of each API call.
    cl_int status;

    // Get the first platform.
    cl_platform_id platform;
    status = clGetPlatformIDs( 1, &platform, &platform_count );
    printf( "%d OpenCL platforms available.\n", (int)platform_count );
    if( status != CL_SUCCESS ) fprintf( stderr, "Error getting platform IDs!\n" );

    // Get the first device.
    cl_device_id device;
    status = clGetDeviceIDs( platform, CL_DEVICE_TYPE_ALL, 1, &device, NULL );
    if( status != CL_SUCCESS ) fprintf( stderr, "Error getting device IDs!\n" );

    // Create a context and associate it with the device.
    cl_context context = clCreateContext( NULL, 1, &device, NULL, NULL, &status );
    if( status != CL_SUCCESS ) fprintf( stderr, "Error creating context!\n" );

    // Create a command queue and associate it with the device.
    cl_command_queue cmdQueue = clCreateCommandQueueWithProperties( context, device, NULL, &status );
    if( status != CL_SUCCESS ) fprintf( stderr, "Error creating the command queue!\n" );

    // Allocate two input buffers and one output buffer.
    cl_mem bufA = clCreateBuffer( context, CL_MEM_READ_WRITE,  datasize * size, NULL, &status );
    cl_mem bufB = clCreateBuffer( context, CL_MEM_READ_WRITE,  datasize, NULL, &status );

    // Transfer data from host arrays into the input buffers.
    status = clEnqueueWriteBuffer( cmdQueue, bufA, CL_FALSE, 0, datasize * size, a, 0, NULL, NULL );
    status = clEnqueueWriteBuffer( cmdQueue, bufB, CL_FALSE, 0, datasize, b, 0, NULL, NULL );


    char *programSource;
    size_t program_size;
    FILE *kernel_file = fopen("elimination_kernel.cl", "rb");
    if (!kernel_file) {
        printf("Failed to load kernel\n");
        return 1;
    }

    fseek(kernel_file, 0, SEEK_END);
    program_size = ftell(kernel_file);
    rewind(kernel_file);
    programSource = (char*)malloc(program_size + 1);
    programSource[program_size] = '\0';
    fread(programSource, sizeof(char), program_size, kernel_file);
    fclose(kernel_file);

    // Create a program with source code.
    cl_program program =
        clCreateProgramWithSource( context, 1, (const char **)&programSource, NULL, &status );
    if( status != CL_SUCCESS ) fprintf( stderr, "Error creating program object!\n" );

    // Build the program for the device.
    status = clBuildProgram( program, 1, &device, NULL, NULL, NULL );
    if( status != CL_SUCCESS ) {
        int log_len;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, NULL, NULL, &log_len);
        char *log = (char*)calloc(log_len, sizeof(char));;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_len, log, NULL);
        switch( status ) {
        case CL_COMPILER_NOT_AVAILABLE:
            fprintf( stderr, "Error building program! Compiler not available.\n" );
            break;
        case CL_BUILD_PROGRAM_FAILURE:
            fprintf( stderr, "Error building program! Build failure.\n" );
            break;
        default:
            fprintf( stderr, "Error building program! *unknown reason*\n" );
            break;
        }
        fprintf( stderr, "%s\n", log );
    }

    // Create the chunk elimination kernel.
    cl_kernel kernel = clCreateKernel( program, "chunk_elimination", &status );
    if( status != CL_SUCCESS ) fprintf( stderr, "Error creating kernel!\n" );

    // Set the kernel arguments.
    status = clSetKernelArg( kernel, 0, sizeof(cl_mem), &bufA );
    status = clSetKernelArg( kernel, 1, sizeof(cl_mem), &bufB );
    status = clSetKernelArg( kernel, 2, sizeof(cl_int), &size );

    // Define an index space of work items for execution.
    size_t indexSpaceSize[1], workGroupSize[1];
    indexSpaceSize[0] = size;
    workGroupSize[0] = 10;

    // Execute the kernel.
    status =
        clEnqueueNDRangeKernel(
          cmdQueue, kernel, 2, NULL, indexSpaceSize, NULL, 0, NULL, NULL);
    if( status != CL_SUCCESS ) fprintf( stderr, "Error queuing kernel execution! Error code: %d\n", status );

    // Read the device output buffers to the host output arrays.
    status = clEnqueueReadBuffer( cmdQueue, bufA, CL_TRUE, 0, datasize, a, 0, NULL, NULL );
    status = clEnqueueReadBuffer( cmdQueue, bufB, CL_TRUE, 0, datasize, b, 0, NULL, NULL );

    enum GaussianResult return_code = gaussian_success; // back_substitution( size, a, b );

    // Free OpenCL resources.
    clReleaseKernel( kernel );
    clReleaseProgram( program );
    clReleaseCommandQueue( cmdQueue );
    clReleaseMemObject( bufA );
    clReleaseMemObject( bufB );
    clReleaseContext( context );

    return return_code;
}

int main( int argc, char *argv[] )
{
    FILE   *input_file;
    size_t  size;

    if( argc < 2 ) {
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

    // Do the calculations.
    Timer stopwatch;
    Timer_initialize( &stopwatch );
    Timer_start( &stopwatch );
    enum GaussianResult result = gaussian_solve( size, a, b );
    Timer_stop( &stopwatch );

    // Display the results.
    switch( result ) {
    case gaussian_success:
        printf( "\nSolution is\n" );
        for( size_t i = 0; i < size; ++i ) {
            printf( " x[%4zu] = %9.5f\n", i, b[i] );
        }
        printf( "Execution time = %ld milliseconds\n", Timer_time( &stopwatch ) );
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
