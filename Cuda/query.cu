
/*
 * Based on a program from Chapter 3 of "CUDA by Example" by Jason Sanders and Edward Kandrot.
 * (C) 2011 by NVIDIA. ISBN=0-13-138768-5.
 */

#include <stdio.h>
#include <stdlib.h>
#include "cuda.h"
#include "driver_types.h"

static void display_boolean( int value )
{
    if( value )
        printf( "TRUE" );
    else
        printf( "FALSE" );
    printf( "\n" );
}


int main( void )
{
    cudaDeviceProp prop;
    int count;
    cudaError_t error_code;

    // cuInit( 0 );

    error_code = cudaGetDeviceCount( &count );
    switch( error_code ) {
    case cudaErrorInsufficientDriver:
    	 printf( "Unable to query devices; driver unavailable!\n" );
	 return EXIT_FAILURE;

    case cudaErrorNoDevice:
    	 printf( "No devices found!\n" );
	 return EXIT_FAILURE;

    case cudaSuccess:
    	 printf( "cudaGetDeviceCount succeeded, %d devices found!\n", count );
	 break;
    }

    for( int i = 0; i < count; ++i ) {
        cudaGetDeviceProperties( &prop, i );
        printf( "\n" );
        printf( "========\n" );
        printf( "Device %d\n", i );
        printf( "========\n" );
        printf( "--- General Information\n" );
        printf( "\tName: %s\n", prop.name );
        printf( "\tIntegrated GPU: " ); display_boolean( prop.integrated );
        printf( "\tCompute capability: %d.%d\n", prop.major, prop.minor );
        printf( "\tClock rate: %d\n", prop.clockRate );
        printf( "\tDevice copy overlap: " ); display_boolean( prop.deviceOverlap );
        printf( "\tKernel execution timeout: " ); display_boolean( prop.kernelExecTimeoutEnabled );

        printf( "\n--- Memory Information\n" );
        printf( "\tTotal global memory: %ld bytes\n", prop.totalGlobalMem );
        printf( "\tTotal constant memory: %ld bytes\n", prop.totalConstMem );
        printf( "\tMaximum memory pitch: %ld bytes\n", prop.memPitch );
        printf( "\tTexture alignment: %ld\n", prop.textureAlignment );
        printf( "\tCan map host memory: " ); display_boolean( prop.canMapHostMemory );

        printf( "\n--- Multiprocessor Information\n" );
        printf( "\tMultiprocessor count: %d\n", prop.multiProcessorCount );
        printf( "\tConcurrent kernels: " ); display_boolean( prop.concurrentKernels );
        printf( "\tShared memory per block: %ld bytes\n", prop.sharedMemPerBlock );
        printf( "\tRegisters per block: %d\n", prop.regsPerBlock );
        printf( "\tThreads in warp: %d\n", prop.warpSize );
        printf( "\tMaximum threads per block: %d\n", prop.maxThreadsPerBlock );
        printf( "\tMaximum block dimensions: (%d, %d, %d)\n",
                prop.maxThreadsDim[0],
                prop.maxThreadsDim[1],
                prop.maxThreadsDim[2] );
        printf( "\tMaximum grid dimensions: (%d, %d, %d)\n",
                prop.maxGridSize[0],
                prop.maxGridSize[1],
                prop.maxGridSize[2] );
        printf( "\n" );
    }
    return EXIT_SUCCESS;
}
