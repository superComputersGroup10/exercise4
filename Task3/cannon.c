#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main (int argc, char **argv) {
    FILE *fp;
    double **A = NULL, **B = NULL, **C = NULL, *C_array = NULL;
    double *A_local_block = NULL, *B_local_block = NULL, *C_local_block = NULL;
    int A_rows, A_columns, A_local_block_rows, A_local_block_columns, A_local_block_size;
    int B_rows, B_columns, B_local_block_rows, B_local_block_columns, B_local_block_size;
    int rank, size, sqrt_size, matrices_a_b_dimensions[4];
    MPI_Comm cartesian_grid_communicator, row_communicator, column_communicator;
    MPI_Status status;

    // used to manage the cartesian grid
    int dimensions[2], periods[2], coordinates[2], remain_dims[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start_total = MPI_Wtime();
    /* For square mesh */
    sqrt_size = (int)sqrt((double) size);
    if(sqrt_size * sqrt_size != size){
        if( rank == 0 ) perror("need to run mpiexec with a perfect square number of processes\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // create a 2D cartesian grid
    dimensions[0] = dimensions[1] = sqrt_size;
    periods[0] = periods[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &cartesian_grid_communicator);
    MPI_Cart_coords(cartesian_grid_communicator, rank, 2, coordinates);

    // create a row communicator
    remain_dims[0] = 0; //ignore the 0th direction
    remain_dims[1] = 1; //take all processes in 1st direction
    MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &row_communicator);

    // create a column communicator
    remain_dims[0] = 1; //take all processes in 0th direction
    remain_dims[1] = 0; //ignore the 1st direction
    MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &column_communicator);


    // READ MATRICES
    int distribs[2], dargs[2], psizes[2];
    MPI_Datatype filetype_A, filetype_B;
    MPI_Offset disp;
    MPI_File fh;

    distribs[0] = MPI_DISTRIBUTE_BLOCK; /* block distribution */
    distribs[1] = MPI_DISTRIBUTE_BLOCK; /* block distribution */
    dargs[0] = MPI_DISTRIBUTE_DFLT_DARG; /* default block size */
    dargs[1] = MPI_DISTRIBUTE_DFLT_DARG; /* default block size */
    psizes[0] = sqrt_size; /* no. of processes in vertical dimension of process grid */
    psizes[1] = sqrt_size; /* no. of processes in horizontal dimension of process grid */

    disp = 2 * sizeof(int); // First two integers indicate #rows and #columns

    // Read Matrix A
    MPI_File_open(MPI_COMM_WORLD, argv[1],  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_all(fh, &matrices_a_b_dimensions[0], 2, MPI_INT, MPI_STATUS_IGNORE);

    MPI_Type_create_darray(size, rank, 2, &matrices_a_b_dimensions[0], distribs, dargs,
            psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype_A);
    MPI_Type_commit(&filetype_A);
    MPI_File_set_view(fh, disp, MPI_DOUBLE, filetype_A, "native", MPI_INFO_NULL);

    A_rows = matrices_a_b_dimensions[0];
    A_columns = matrices_a_b_dimensions[1];
    A_local_block_rows = (int) (A_rows / sqrt_size);
    A_local_block_columns = (int) (A_columns / sqrt_size);
    A_local_block_size = A_local_block_rows * A_local_block_columns;
    A_local_block = (double *) malloc (A_local_block_size * sizeof(double));
    MPI_File_read_all(fh, A_local_block, A_local_block_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    // Read Matrix B
    MPI_File_open(MPI_COMM_WORLD, argv[2],  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_all(fh, &matrices_a_b_dimensions[2], 2, MPI_INT, MPI_STATUS_IGNORE);

    MPI_Type_create_darray(size, rank, 2, &matrices_a_b_dimensions[2], distribs, dargs,
            psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype_B);
    MPI_Type_commit(&filetype_B);
    MPI_File_set_view(fh, disp, MPI_DOUBLE, filetype_B, "native", MPI_INFO_NULL);

    B_rows = matrices_a_b_dimensions[2];
    B_columns = matrices_a_b_dimensions[3];
    B_local_block_rows = (int) B_rows / sqrt_size;
    B_local_block_columns = (int) B_columns / sqrt_size;
    B_local_block_size = B_local_block_rows * B_local_block_columns;
    B_local_block = (double *) malloc (B_local_block_size * sizeof(double));
    MPI_File_read_all(fh, B_local_block, B_local_block_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    // need to check that the multiplication is possible given dimensions
    // matrices_a_b_dimensions[0] = row size of A
    // matrices_a_b_dimensions[1] = column size of A
    // matrices_a_b_dimensions[2] = row size of B
    // matrices_a_b_dimensions[3] = column size of B
    if(matrices_a_b_dimensions[1] != matrices_a_b_dimensions[2]){
        if(rank == 0) fprintf(stderr, "A's column size (%d) must match B's row size (%d)\n",
                matrices_a_b_dimensions[1], matrices_a_b_dimensions[2]);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // this implementation is limited to cases where the matrices can be partitioned perfectly
    if( matrices_a_b_dimensions[0] % sqrt_size != 0
            || matrices_a_b_dimensions[1] % sqrt_size != 0
            || matrices_a_b_dimensions[2] % sqrt_size != 0
            || matrices_a_b_dimensions[3] % sqrt_size != 0 ){
        if(rank == 0) fprintf(stderr, "cannot distribute work evenly among processe\n"
                "all dimensions (A: r:%d c:%d; B: r:%d c:%d) need to be divisible by %d\n",
                matrices_a_b_dimensions[0],matrices_a_b_dimensions[1],
                matrices_a_b_dimensions[2],matrices_a_b_dimensions[3], sqrt_size );
        MPI_Abort(MPI_COMM_WORLD, -1);
    }


    // local metadata for C
    C_local_block = (double *) malloc (A_local_block_rows * B_local_block_columns * sizeof(double));
    // C needs to be initialized at 0 (accumulates partial dot-products)
    int i;
    for(i=0; i < A_local_block_rows * B_local_block_columns; i++){
        C_local_block[i] = 0;
    }



    //Initial matrix alignment for A and B
    int shift_source, shift_destination;
    MPI_Cart_shift(row_communicator, 0, -coordinates[0], &shift_source, &shift_destination);
    MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE,
            shift_destination, 0,
            shift_source, 0, row_communicator, &status);

    MPI_Cart_shift(column_communicator, 0, -coordinates[1], &shift_source, &shift_destination);
    MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE,
            shift_destination, 0,
            shift_source, 0, column_communicator, &status);

    // cannon's algorithm
    int cannon_block_cycle;
    double compute_time = 0, mpi_time = 0, start;   //all ranks declare compute_time, mpi_time
    int C_index, A_row, A_column, B_column;
    for(cannon_block_cycle = 0; cannon_block_cycle < sqrt_size; cannon_block_cycle++){
        // compute partial result for this block cycle
        start = MPI_Wtime();
        for(C_index = 0, A_row = 0; A_row < A_local_block_rows; A_row++){
            for(B_column = 0; B_column < B_local_block_columns; B_column++, C_index++){
                for(A_column = 0; A_column < A_local_block_columns; A_column++){
                    C_local_block[C_index] += A_local_block[A_row * A_local_block_columns + A_column] *
                        B_local_block[A_column * B_local_block_columns + B_column];
                }
            }
        }

        compute_time += MPI_Wtime() - start;        //each rank accumulates the compute_time
        start = MPI_Wtime();
        // rotate blocks horizontally
        MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE,
                (coordinates[1] + sqrt_size - 1) % sqrt_size, 0,
                (coordinates[1] + 1) % sqrt_size, 0, row_communicator, &status);
        // rotate blocks vertically
        MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE,
                (coordinates[0] + sqrt_size - 1) % sqrt_size, 0,
                (coordinates[0] + 1) % sqrt_size, 0, column_communicator, &status);
        mpi_time += MPI_Wtime() - start;
    }


    // Write output matrix
    MPI_Datatype filetype_C;
    int size_C[2];
    size_C[0] = A_rows;
    size_C[1] = B_columns;
    const int C_local_block_size = A_local_block_rows * B_local_block_columns;

    MPI_File_open(MPI_COMM_WORLD, argv[3],  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    // Only rank 0 writes the header
    if(rank == 0) MPI_File_write(fh, size_C, 2, MPI_INT, MPI_STATUS_IGNORE);
    MPI_Type_create_darray(size, rank, 2, size_C, distribs, dargs,
            psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype_C);
    MPI_Type_commit(&filetype_C);
    MPI_File_set_view(fh, disp, MPI_DOUBLE, filetype_C, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, C_local_block, C_local_block_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    double end_total = MPI_Wtime() - start_total;

    // generating output at rank 0
    if (rank == 0) {
        printf("(%d,%d)x(%d,%d)=(%d,%d)\n", A_rows, A_columns, B_rows, B_columns, A_rows, B_columns);
        printf("Computation time: %lf\n", compute_time);
        printf("MPI time:         %lf\n", mpi_time);
        printf("Total time:       %lf\n", end_total);
    }

    if (argc == 5){
        // full array only needed at root
        if(rank == 0){
            C_array = (double *) malloc(sizeof(double) * A_rows * B_columns);
            // allocate output matrix C
            C = (double **) malloc(A_rows * sizeof(double *));
            for(i=0; i<A_rows ;i++){
                C[i] = (double *) malloc(B_columns * sizeof(double));
            }
        }
        // get C parts from other processes at rank 0
        MPI_Gather(C_local_block, A_local_block_rows*B_local_block_columns, MPI_DOUBLE,
                C_array, A_local_block_rows*B_local_block_columns, MPI_DOUBLE,
                0, cartesian_grid_communicator);

        if (rank == 0){
            // convert the ID array into the actual C matrix
            int i, j, k, row, column;
            for (i = 0; i < sqrt_size; i++){  // block row index
                for (j = 0; j < sqrt_size; j++){ // block column index
                    for (row = 0; row < A_local_block_rows; row++){
                        for (column = 0; column < B_local_block_columns; column++){
                            C[i * A_local_block_rows + row] [j * B_local_block_columns + column] =
                                C_array[((i * sqrt_size + j) * A_local_block_rows * B_local_block_columns)
                                + (row * B_local_block_columns) + column];
                        }
                    }
                }
            }


            printf("\nPerforming serial consistency check. Be patient...\n");
            fflush(stdout);
            int pass = 1;
            double temp;

            double * A_array = (double *) malloc(A_rows * A_columns * sizeof(double));
            double * B_array = (double *) malloc(B_rows * B_columns * sizeof(double));

            //  Read again from file
            MPI_File_open(MPI_COMM_SELF, argv[1],  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            MPI_File_read_at(fh, disp, A_array, A_rows * A_columns, MPI_DOUBLE, MPI_STATUS_IGNORE);
            MPI_File_close(&fh);

            MPI_File_open(MPI_COMM_SELF, argv[2],  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            MPI_File_read_at(fh, disp, B_array, B_rows * B_columns, MPI_DOUBLE, MPI_STATUS_IGNORE);
            MPI_File_close(&fh);

            // Checking
            for(i=0; i<A_rows; i++){
                for(j=0; j<B_columns; j++){
                    temp = 0;
                    for(k=0; k<B_rows; k++){
                        temp += A_array[i*A_columns + k] * B_array[k*B_columns + j];
                    }
                    if(temp != C[i][j]){
                        pass = 0;
                    }
                }
            }
            if (pass) printf("Consistency check: PASS\n");
            else printf("Consistency check: FAIL\n");

            // Free memory
            for(int i = 0; i < A_rows; i++)
                free(C[i]);
            free(C);
            free(C_array);
            free(A_array);
            free(B_array);
        }
    }


    // free all memory
    free(A_local_block);
    free(B_local_block);
    free(C_local_block);

    // finalize MPI
    MPI_Finalize();
    return 0;
}
