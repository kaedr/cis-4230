
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#define N 100

__global__ void add( float (*a)[N], float (*b)[N], float (*c)[N] )
{
    int row = blockIdx.x;
    int col = blockIdx.y;
    if( row < N && col < N )
        c[row][col] = a[row][col] + b[row][col];
}


float a[N][N];
float b[N][N];
float c[N][N];

void initialize( )
{
    for( int row = 0; row < N; ++row ) {
        for( int col = 0; col < N; ++col ) {
            a[row][col] = 1.0;
            b[row][col] = 2.0;
        }
    }
}


void check_result( )
{
    for( int row = 0; row < N; ++row ) {
        for( int col = 0; col < N; ++col ) {
            if( fabs( c[row][col] - 3.0 ) > 1E-06 ) {
                printf( "FAILED!\n" );
                goto done;
            }
        }
    }
    printf( "PASSED!\n" );
done: ;
}


int main( void )
{
    float (*dev_a)[N], (*dev_b)[N], (*dev_c)[N];
    cudaError_t result;

    initialize( );

    cudaMalloc( (void **)&dev_a, N * N * sizeof(float) );
    cudaMalloc( (void **)&dev_b, N * N * sizeof(float) );
    cudaMalloc( (void **)&dev_c, N * N * sizeof(float) );

    cudaMemcpy( dev_a, a, N * N * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, N * N * sizeof(float), cudaMemcpyHostToDevice );
    dim3 grid(N, N);
    add<<<grid, 1>>>( dev_a, dev_b, dev_c );
    result = cudaGetLastError( );
    if( result != cudaSuccess ) {
        printf( "Failure invoking the add kernel! (%s)\n", cudaGetErrorString( result ) );
    }
    cudaMemcpy( c, dev_c, N * N * sizeof(float), cudaMemcpyDeviceToHost );

    cudaFree( dev_a );
    cudaFree( dev_b );
    cudaFree( dev_c );

    check_result( );

    return EXIT_SUCCESS;
}
