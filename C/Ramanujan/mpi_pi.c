
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gmp.h>
#include <mpfr.h>
#include <omp.h>
#include <mpi.h>

#define DECIMAL_PRECISION    10000   // Digits of precision

// A struct for memoization of previous calculations
struct Memo {
    // Below are expressed the calculations of these terms
    // for any future k (which we'll call n)
    unsigned int k;
    // (4n)! = (4k!) * 4n * 4(n-1) ... 4(k+1)
    mpz_t four_k_factorial;
    // n!^4 = k!^4 * n^4 * (n-1)^4 ... (k+1)^4
    mpz_t k_factorial_fourth;
    // 26390n = 26390k * 26390(n-k)
    mpz_t k_26390;
    // 396^4n = 396^4k * 396^4(n-k)
    mpz_t three_ninety_six_to_four_k;
    // Store the sqrt(8) /9801 factor
    mpfr_t factor;
    // Try to avoid reinitializing the term each loop
    mpfr_t term;
};

// Intialize the Memo to a useful starting state
void memo_init(struct Memo *k_memo, unsigned int desired_precision) {
    k_memo->k = 0;
    // Init the oversized ints
    mpz_init_set_ui( k_memo->four_k_factorial, 1UL );
    mpz_init_set_ui( k_memo->k_factorial_fourth, 1UL );
    mpz_init_set_ui( k_memo->k_26390, 0UL );
    mpz_init_set_ui( k_memo->three_ninety_six_to_four_k, 1UL );

    // Init the giant floats
    mpfr_init2( k_memo->factor, desired_precision );
    mpfr_init2( k_memo->term, desired_precision );
    // And populate initial value of factor
    mpfr_sqrt_ui( k_memo->factor, 8UL, MPFR_RNDN );
    mpfr_div_ui( k_memo->factor, k_memo->factor, 9801UL, MPFR_RNDN );
}

void memo_destroy(struct Memo *k_memo) {
    // Clean up oversized ints
    mpz_clear( k_memo->four_k_factorial );
    mpz_clear( k_memo->k_factorial_fourth );
    mpz_clear( k_memo->k_26390 );
    mpz_clear( k_memo->three_ninety_six_to_four_k );

    // Clean up giant floats
    mpfr_clear( k_memo->term );
    mpfr_clear( k_memo->factor );
}

void print_memo(struct Memo *k_memo) {
    printf ("Memo { k : %u || ", k_memo->k );
    gmp_printf ("26390k : %Zd || ", k_memo->k_26390);
    gmp_printf ("396^4k : %Zd || ", k_memo->three_ninety_six_to_four_k);
    gmp_printf ("4k! : %Zd || ", k_memo->four_k_factorial);
    gmp_printf ("k!^4 : %Zd }\n", k_memo->k_factorial_fourth);

}

// (4n)! = (4k!) * 4n * 4n-1 ... 4k+1
void next_four_k_factorial(mpz_t current, unsigned int k, unsigned int n) {
    for (unsigned int i = n; i > k; --i) {
        unsigned int fact_term = 4*i;
        mpz_mul_ui(current, current, fact_term);
        mpz_mul_ui(current, current, fact_term - 1);
        mpz_mul_ui(current, current, fact_term - 2);
        mpz_mul_ui(current, current, fact_term - 3);
    }
}

 // n!^4 = k!^4 * n^4 * (n-1)^4 ... (k+1)^4
void next_k_factorial_fourth(mpz_t current, unsigned int k, unsigned int n) {
    mpz_t n_minus_x_fourth;
    mpz_init( n_minus_x_fourth );
    for (unsigned int i = n; i > k; --i) {
        mpz_ui_pow_ui(n_minus_x_fourth, i, 4UL); // (n-x)^4
        mpz_mul(current, current, n_minus_x_fourth);
    }

    // Free up temp
    mpz_clear( n_minus_x_fourth );
}

// 26390n = 26390k + 26390(n-k)
void next_k_26390(mpz_t current, unsigned int k, unsigned int n) {
    mpz_add_ui(current, current, 26390 * (n-k));
}

// 396^4n = 396^4k * 396^4(n-k)
void next_three_ninety_six_to_four_k(mpz_t current, unsigned int k, unsigned int n) {
    mpz_t three_ninety_six_to_four_n_minus_k;
    mpz_init( three_ninety_six_to_four_n_minus_k );
    mpz_ui_pow_ui(three_ninety_six_to_four_n_minus_k, 396UL, 4UL*(n-k));
    mpz_mul(current, current, three_ninety_six_to_four_n_minus_k);

    // Free up temp
    mpz_clear( three_ninety_six_to_four_n_minus_k );
}

void advance_to_n(struct Memo *k_memo, unsigned int n) {
    next_four_k_factorial(k_memo->four_k_factorial, k_memo->k, n);
    next_k_factorial_fourth(k_memo->k_factorial_fourth, k_memo->k, n);
    next_k_26390(k_memo->k_26390, k_memo->k, n);
    next_three_ninety_six_to_four_k(k_memo->three_ninety_six_to_four_k, k_memo->k, n);
    k_memo->k = n;
}

void update_accumulator(struct Memo *k_memo, mpfr_t accumulator) {
    // Because our first action is to set the value of the term,
    // We avoid having to do any other cleanup each loop
    // Set it to 26390k + 1103
    mpfr_set_z( k_memo->term, k_memo->k_26390, MPFR_RNDN );
    mpfr_add_ui( k_memo->term, k_memo->term, 1103UL, MPFR_RNDN );
    // Divide by 396^4k
    mpfr_div_z( k_memo->term, k_memo->term, k_memo->three_ninety_six_to_four_k, MPFR_RNDN );
    // Multiply by 4k!
    mpfr_mul_z( k_memo->term, k_memo->term, k_memo->four_k_factorial, MPFR_RNDN );
    // Divide by k!^4
    mpfr_div_z( k_memo->term, k_memo->term, k_memo->k_factorial_fourth, MPFR_RNDN );


    // multiply the factor into our term
    mpfr_mul( k_memo->term, k_memo->term, k_memo->factor, MPFR_RNDN );

    // Update accumulator
    mpfr_add( accumulator, accumulator, k_memo->term, MPFR_RNDN );
}

int check_results(mpfr_t accumulator, int precision) {

    FILE * pi_file = fopen("pi.txt", "r");
    if (pi_file == NULL) {
        printf("can't find pi.txt");
        return EXIT_FAILURE;
    }
    char * buffer;
    char * calc_buffer = NULL;
    long   numbytes;
    // Get the number of bytes
    fseek(pi_file, 0L, SEEK_END);
    numbytes = ftell(pi_file);

    // reset the file position indicator to
    // the beginning of the file
    fseek(pi_file, 0L, SEEK_SET);

    // grab sufficient memory for the
    // buffer to hold the text
    buffer = (char*)calloc(numbytes, sizeof(char));

    // memory error
    if(buffer == NULL)
        return EXIT_FAILURE;

    // printf( "reading file\n" );
    // copy all the text into the buffer
    fread( buffer, sizeof(char), numbytes, pi_file );
    fclose( pi_file );

    // printf( "Getting pi string\n" );
    // Write our calulated pi to the buffer
    long i = 1;
    long * i_ptr = &i;
    calc_buffer = mpfr_get_str( calc_buffer, i_ptr, 10, precision, accumulator, MPFR_RNDN );


    // Check our accuracy
    for (unsigned int overlap = 0; overlap <= precision; ++overlap) {
        // printf("Comparing %c to %c\n", buffer[overlap], calc_buffer[overlap]);
        if ( buffer[overlap] != calc_buffer[overlap] ) {
            printf("Comparing %c to %c\n", buffer[overlap], calc_buffer[overlap]);
            printf("Calculation accurate to %u decimal places\n", overlap);
            break;
        }
    }

    // free the memory we used for the buffer
    free(buffer);
    mpfr_free_str(calc_buffer);

    return EXIT_SUCCESS;
}

void thread_work(
        mpfr_t accumulators[], struct Memo k_memos[], unsigned int threads,
        unsigned int iterations, int number_of_nodes, int local_rank
    ) {
    #pragma omp parallel
    {

        double start_time;
        double checkpoint;
        unsigned int TID = omp_get_thread_num();
        // Make sure we start in the right place
        unsigned int thread_local_offset = (local_rank * threads) + TID;
        // Increment by threads * nodes
        unsigned int increment = number_of_nodes * threads;
        // How often to output
        unsigned int interval = 4000;

        if ( local_rank == 0 && TID == 0 ) {
            // Replacing clock() with omp_get_wtime() to correctly handle threading
            start_time = omp_get_wtime();
        }

        for (unsigned int k = thread_local_offset; k <= iterations; k += increment) {
            if ( local_rank == 0 && TID == 0 && k > 0 && k % interval == 0) {
                checkpoint = omp_get_wtime() - start_time;
                // Because clock counts cpu time, it advances
                printf("%u iterations complete in %fs\n", k, checkpoint);
            }
            // printf("Thread %u executing iteration %u\n", TID, k);
            advance_to_n( &k_memos[TID], k );
            // print_memo(&k_memo);
            update_accumulator( &k_memos[TID], accumulators[TID]);
        }
        if ( local_rank == 0 && TID == 0 ) {
            checkpoint = omp_get_wtime() - start_time;
            // Because clock counts cpu time, it advances
            printf("%u iterations complete in %fs\n", iterations, checkpoint);
        }
    }
}

int main( int argc, char *argv[] ) {
    int number_of_nodes;
    int local_rank;
    int source_node;
    const int destination_node = 0;
    const int tag_value = 8675309;
    MPI_Status status;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_nodes );
    MPI_Comm_rank( MPI_COMM_WORLD, &local_rank );

    // Handle command line arg for changing precision
    char *char_pointer;
    unsigned long precision = DECIMAL_PRECISION;
    if (argc > 1) {
        errno = 0;
        long converted = strtol(argv[1], &char_pointer, 10);

        if (errno == 0 && *char_pointer == '\0' && converted > 0 && converted <= __LONG_MAX__) {
            precision = converted;
        }
    }
    // Handle command line arg for changing thread pool
    if (argc > 2) {
        errno = 0;
        long converted = strtol(argv[2], &char_pointer, 10);

        if (errno == 0 && *char_pointer == '\0' && converted > 0 && converted <= __LONG_MAX__) {
            omp_set_num_threads(converted);
        }
    }

    // How many bits are needed to achieve the desired decimal precision?
    unsigned int desired_precision = (int)ceil( precision * ( log( 10 ) / log( 2 ) ) );

    // Find out how many threads we'll be working with
    unsigned int threads = omp_get_max_threads();

    if ( local_rank == 0 && number_of_nodes > threads ) {
        printf("You'll need to do something different with accumulators for this to work");
        MPI_Finalize( );
        return EXIT_FAILURE;
    }

    // Set up our accumulators and memos
    mpfr_t accumulators[threads];
    struct Memo k_memos[threads];

    for (int i = 0; i < threads; ++i) {
        mpfr_init2( accumulators[i], desired_precision );
        mpfr_set_d( accumulators[i], 0.0, MPFR_RNDN );
        memo_init( &k_memos[i], desired_precision );
    }

    // For practical purposes we seem to get 7 digits precision per iteration
    // The +1 is to ensure we don't get shortchanged by the integer division
    unsigned int iterations = (precision / 7) + 1;

    if ( local_rank == 0 ) {
        printf("Making %u iterations across %u threads\n", iterations, threads);
    }

    // Do the work
    thread_work(accumulators, k_memos, threads, iterations, number_of_nodes, local_rank);

    // collect thread work
    for (int i = 1; i < threads; ++i) {
        mpfr_add(accumulators[0], accumulators[0], accumulators[i], MPFR_RNDN );
    }

    // Collect node work
    // Worth noting, there's some of this that doesn't strictly need to happen on all nodes
    // But for now I'm not worried about the time wasted by a few operations at collection.
    // 1. Prep for serialization
    char * send_buffer;
    size_t send_buffer_size;
    FILE * mpfr_t_stream = open_memstream(&send_buffer, &send_buffer_size);
    // 2. Serialize to stream
    mpfr_fpif_export(mpfr_t_stream, accumulators[0]);
    fclose(mpfr_t_stream);

    // 3. Prep a space to gather info about the length of what we'll receive
    size_t buffer_sizes[number_of_nodes];
    // 4. Gather up the length info that'll need to receive the accumulators
    MPI_Gather(&send_buffer_size, 1, MPI_UNSIGNED_LONG, buffer_sizes, 1, MPI_UNSIGNED_LONG, destination_node, MPI_COMM_WORLD);

    int result;
    // I went with != here so that the sends come before receives in code, for my own mental clarity
    if ( local_rank != 0 ) {
        MPI_Send(send_buffer, send_buffer_size, MPI_BYTE, destination_node, tag_value, MPI_COMM_WORLD);
    } else {
        // Main node collects
        MPI_Status status;
        for ( source_node = 1; source_node < number_of_nodes; ++source_node ) {
            char byte_bin[buffer_sizes[source_node]];
            MPI_Recv(byte_bin, buffer_sizes[source_node], MPI_BYTE, source_node, tag_value, MPI_COMM_WORLD, &status);
            mpfr_t_stream = fmemopen(byte_bin, buffer_sizes[source_node], "rb");
            // I'm reusing the accumulators array here, which is safe because of an earlier check
            mpfr_fpif_import(accumulators[source_node], mpfr_t_stream);
            mpfr_add(accumulators[0], accumulators[0], accumulators[source_node], MPFR_RNDN );
        }

        // account for the 1/pi thing
        mpfr_ui_div( accumulators[0], 1UL, accumulators[0], MPFR_RNDN );

        // Print the answer.
        printf( "pi = " );
        mpfr_out_str( stdout, 10, 50, accumulators[0], MPFR_RNDN );
        printf( "\n" );

        result = check_results(accumulators[0], precision);
    }



    for (int i = 0; i < threads; ++i) {
        // Release resources associated with the multi-precision floats.
        mpfr_clear( accumulators[i] );
        // Clear the memo struct
        memo_destroy(&k_memos[i]);
    }

    MPI_Finalize( );
    return result;
}
