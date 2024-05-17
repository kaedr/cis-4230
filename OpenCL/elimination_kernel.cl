#define MATRIX_GET( matrix, size, row, column )        ( (matrix)[(row)*(size) + (column)] )
#define MATRIX_PUT( matrix, size, row, column, value ) ( (matrix)[(row)*(size) + (column)] = (value) )

__kernel void chunk_elimination(__global double *a, __global double *b, unsigned int size, unsigned int i) {
    int my_id = get_global_id( 0 );
    size_t         k;
    double  m;
    int j = i + 1 + my_id;

    m = MATRIX_GET( a, size, j, i ) / MATRIX_GET( a, size, i, i );
    for( k = 0; k < size; ++k )
        MATRIX_PUT( a, size, j, k, MATRIX_GET( a, size, j, k ) - m * MATRIX_GET( a, size, i, k ) );
    b[j] -= m * b[i];
}
