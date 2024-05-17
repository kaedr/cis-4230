#define MATRIX_GET( matrix, size, row, column )        ( (matrix)[(row)*(size) + (column)] )
#define MATRIX_PUT( matrix, size, row, column, value ) ( (matrix)[(row)*(size) + (column)] = (value) )

__kernel void chunk_elimination(__global double *a, __global double *b, unsigned int size) {

    unsigned int    i, j, k, l;
    double          temp, m;
    int my_id = get_global_id( 0 );

    for( i = 0; i < size - 1; ++i ) {
        barrier(CLK_GLOBAL_MEM_FENCE);
        if (my_id > size - i ) continue;
        if (my_id == 0) {
            // Find the row with the largest value of |a[j][i]|, j = i, ..., n - 1
            k = i;
            m = fabs( MATRIX_GET(a, size, i, i) );
            for( j = i + 1; j < size; ++j ) {
                if( fabs( MATRIX_GET(a, size, j, i) ) > m ) {
                    k = j;
                    m = fabs( MATRIX_GET(a, size, j, i) );
                }
            }

            // Check for |a[k][i]| zero.
            // TODO: The value 1.0E-6 is arbitrary. A more disciplined value should be used.
            if( fabs( MATRIX_GET(a, size, k, i) ) <= 1.0E-6 ) {
                //free( temp_array );
                return;
            }

            // Exchange row i and row k, if necessary.
            if( k != i ) {
                // memcpy( temp_array, a[i], size * sizeof( double ) );
                // memcpy( a[i], a[k], size * sizeof( double ) );
                // memcpy( a[k], temp_array, size * sizeof( double ) );
                for ( l = 0; l < size; l++) {
                    temp = MATRIX_GET(a, size, i, l);
                    MATRIX_PUT(a, size, i, l, MATRIX_GET(a, size, k, l));
                    MATRIX_PUT(a, size, k, l, temp);
                }

                // Exchange corresponding elements of b.
                temp = b[i];
                b[i] = b[k];
                b[k] = temp;
            }
            continue;
        }

        barrier(CLK_GLOBAL_MEM_FENCE);

        j = i + 1 + my_id;
        m = MATRIX_GET( a, size, j, i ) / MATRIX_GET( a, size, i, i );
        for( k = 0; k < size; ++k )
            MATRIX_PUT( a, size, j, k, MATRIX_GET( a, size, j, k ) - m * MATRIX_GET( a, size, i, k ) );
        b[j] -= m * b[i];
    }
}
