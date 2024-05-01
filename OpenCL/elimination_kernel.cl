__kernel void chunk_elimination(__global double *a, __global double *b, int size) {
    int idx = get_global_id( 0 );
    int offset = idx * size;
    for (int i = 0; i < size; ++i) {
        b[idx] += a[offset + i];
    }
    // b[idx] = (double)idx;
    /* code */
}
