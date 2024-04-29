
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

__global__ void add( int a, int b, int *c )
{
    *c = a + b;
}

int main( void )
{
    int  c;
    int *dev_c;
    cudaError_t result;

    result = cudaMalloc( (void **)&dev_c, sizeof(int) );
    if( result != cudaSuccess ) {
        printf( "Failure to allocate memory on the device!\n" );
    }

    add<<<1, 1>>>( 2, 7, dev_c );
    result = cudaGetLastError( );
    if( result != cudaSuccess ) {
        printf( "Failure invoking the add kernel! (%s)\n", cudaGetErrorString( result ) );
    }

    result = cudaMemcpy( &c, dev_c, sizeof(int), cudaMemcpyDeviceToHost );
    if( result != cudaSuccess ) {
        printf( "Failure to copy result from the device to the host!\n" );
    }

    cudaFree( dev_c );
    printf( "c = %d\n", c );
    return EXIT_SUCCESS;
}
