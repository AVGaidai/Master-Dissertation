/** \file */
/*
 * GAUSS ELIMINATION (MPI-FT VERSION) BY ANATOLY GAIDAI
 * 2018
 */
#include <stdio.h>
#include <stdlib.h>   /* malloc() */

#include <signal.h>

#include <unistd.h>   /* sleep()  */

#include <string.h>   /* memcpy() */

#include <sys/time.h> /* gettimeofday () */

#include <time.h>

#include <math.h>

#include <mpi.h>
#include <mpi-ext.h>


#define CHANCE 0.0005


MPI_Comm COMM, NEWCOMM;
MPI_Errhandler new_eh;

struct timeval start, end, t;

int mpi_mcw_rank, mpi_mcw_size;
int flg = 0;


char * get_str_failed_procs(MPI_Comm comm, MPI_Group f_group)
{
    int f_size, i, c_size;
    MPI_Group c_group;
    int *failed_ranks = NULL;
    int *comm_ranks   = NULL;
    int rank_len = 7;
    char * ranks_failed = NULL;

    MPI_Group_size(f_group, &f_size);

    if( f_size <= 0 ) {
        ranks_failed = strdup("None");
    } else {
        MPI_Comm_group(comm, &c_group);
        MPI_Comm_size( comm, &c_size);

        failed_ranks = (int *)malloc(f_size * sizeof(int));
        comm_ranks   = (int *)malloc(f_size * sizeof(int));
        for( i = 0; i < f_size; ++i) {
            failed_ranks[i] = i;
        }

        MPI_Group_translate_ranks(f_group, f_size, failed_ranks,
                                  c_group, comm_ranks);

        ranks_failed = (char *)malloc(sizeof(char) * (rank_len) * f_size + 1);
        for( i = 0; i < f_size; ++i) {
            /*
            printf("%2d of %2d) Error Handler: %2d / %2d Failed Rank %3d\n",
                   mpi_mcw_rank, mpi_mcw_size, i, f_size, comm_ranks[i]);
            */
            if( i+1 == f_size ) {
                sprintf((ranks_failed+(i*rank_len)), "%c%5d.%c",
                        (0 == i ? ' ' : ','), comm_ranks[i], '\0');
            }
            else if( 0 == i ) {
                sprintf(ranks_failed, "  %5d", comm_ranks[i]);
            }
            else {
                sprintf((ranks_failed+(i*rank_len)), ", %5d", comm_ranks[i]);
            }
        }

        MPI_Group_free(&c_group);
    }

    if( NULL != failed_ranks ) {
        free(failed_ranks);
    }
    if( NULL != comm_ranks ) {
        free(comm_ranks);
    }

    return ranks_failed;
}


void mpi_error_handler(MPI_Comm *comm, int *error_code, ...)
{
    MPI_Group f_group;
    int num_failures;
    int loc_size;
    char * ranks_failed = NULL;

    MPI_Comm_size(*comm, &loc_size);

    if( *error_code == MPI_ERR_PROC_FAILED ) {
        printf("Process %d\n", mpi_mcw_rank);
        /* Access the local list of failures */
        MPIX_Comm_failure_ack(*comm);
        MPIX_Comm_failure_get_acked(*comm, &f_group);

        /* Get the number of failures */
        MPI_Group_size(f_group, &num_failures);

        ranks_failed = get_str_failed_procs(*comm, f_group);

        printf("%2d of %2d) Error Handler: (Comm = %s) %3d Failed Ranks: %s\n",
               mpi_mcw_rank, mpi_mcw_size,
               (mpi_mcw_size == loc_size ? "MCW" : "Subcomm"),
               num_failures, ranks_failed);
        free(ranks_failed);

        MPIX_Comm_agree(COMM, &loc_size);
        MPIX_Comm_revoke(COMM);
        MPIX_Comm_shrink(COMM, &NEWCOMM);
        MPI_Comm_free(&COMM);
        COMM = NEWCOMM;
        MPI_Comm_set_errhandler(COMM, new_eh);

        MPI_Comm_rank(COMM, &mpi_mcw_rank);
    } else {
        return;
    }

    /* Introduce a small delay to aid debugging */
    fflush(NULL);
    
    return;
}



/**
 * \brief Print content of matrix 'matrix'.
 *
 * \param matrix is pointer on matrix.
 * \param SIZE_X is number of rows in the matrix.
 * \param SIZE_Y is number of columns in the matrix.
 *
 * \return no value.
 */
void print_matrix(double *matrix, int SIZE_X, int SIZE_Y)
{
    int i, j;

    printf("Matrix:\n");
    for (i = 0; i < SIZE_X; ++i) {
	for (j = 0; j < SIZE_Y; ++j) {
	    printf("%.2f\t", matrix[i * SIZE_Y + j]);
	}
	printf("\n");
    }
}


/**
 * \brief Partition of the input matrix to all processes.
 *
 * \param fp is file descriptor of the input matrix.
 * \param matrix is pointer to the matirx.
 * \param ROWS is pointer to the number of rows in the matrix.
 * \param COLUMNS is pointer to the number of columns in the matrix.
 *

 * \return no value.
 */
void MPI_Matrix_Partition(FILE *fp, double **matrix,
                          int *ROWS, int *COLUMNS)
{
    int i, j;
    double *ptr;
    
    /* Root process */
    if (mpi_mcw_rank == 0) {
        /* Reading matrix size */
        fscanf(fp, "%d %d", ROWS, COLUMNS);
        
        *matrix = (double *) malloc((*COLUMNS) * (*ROWS) * sizeof(double));
        ptr = *matrix;
        
        for (i = 0; i < *ROWS; ++i) {
            for (j = 0; j < *COLUMNS; ++j) {
                fscanf(fp, "%lf", &ptr[i * (*COLUMNS) + j]);
            }
        }
    }

    MPI_Bcast((void *) COLUMNS, 1, MPI_INT, 0, COMM);
    MPI_Bcast((void *) ROWS, 1, MPI_INT, 0, COMM);

    if (mpi_mcw_rank != 0) {
        *matrix = (double *) malloc((*COLUMNS) * (*ROWS) *sizeof(double));
    }
    
    MPI_Bcast((void *) (*matrix), (*COLUMNS) * (*ROWS),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


int Failure_gen(double chance)
{
    double x;
    
    gettimeofday(&t, NULL);
    srand((unsigned int) t.tv_usec * (mpi_mcw_rank + 1));
    x = rand() / (double) RAND_MAX;
    
    if (x <= chance) {
        printf("%f\n", x);
        return 1;
    }
    return 0;
}



/**
 * \brief Matrix transformation to lower-triangular.
 *
 * \param matrix is pointer to the matrix.
 * \param ROWS is number of rows in the matrix.
 * \param COLUMNS is number of columns in the matrix.
 *
 * \return no value.
 */
void MPI_Gauss_Forward(double *matrix, int ROWS, int COLUMNS)
{
    int curr_row = 0;
    int i, j, l, proc_rows, block_size;
    double k;
    
    srand(time(NULL));
    for ( ; curr_row < ROWS; ++curr_row) {
        k = matrix[curr_row * COLUMNS + curr_row];
        for (i = curr_row + 1; i < COLUMNS; ++i) {
            matrix[curr_row * COLUMNS + i] /= k;
        }
        matrix[curr_row * COLUMNS + curr_row] = 1.00;

        if (curr_row == ROWS - 1) continue;
        
        proc_rows = ROWS - (curr_row + 1);
        block_size = proc_rows / mpi_mcw_size;
        l = curr_row + 1 + block_size * mpi_mcw_rank;
        if (mpi_mcw_rank == mpi_mcw_size - 1) {
            block_size += proc_rows % mpi_mcw_size;
        }
    
        for (i = 0; i < block_size; ++i) {
            k = 0.00 - matrix[(l + i) * COLUMNS + curr_row]; 
            for (j = curr_row + 1; j < COLUMNS; ++j) {
                matrix[(l + i) * COLUMNS + j] +=
                    matrix[curr_row * COLUMNS + j] * k;
            }
            matrix[(l + i) * COLUMNS + curr_row] = 0.00;
        }
        /* sleep(mpi_mcw_rank); */
        /* printf("p%d\n", mpi_mcw_rank); */
        /* print_matrix(matrix, ROWS, COLUMNS); */
        
        if (mpi_mcw_rank == 0) {
            if (Failure_gen(CHANCE)) {
                printf("(Forward) N %d \t STEP %d\n",
                       mpi_mcw_size - 1, curr_row); 
                raise(SIGKILL);
            }
        }
        
        MPI_Barrier(COMM);
        /* printf("message 1\n"); */
        i = mpi_mcw_size;
        MPI_Comm_size(COMM, &mpi_mcw_size);
        if (mpi_mcw_size < i) {
            --curr_row;
            MPI_Comm_rank(COMM, &mpi_mcw_rank);
            continue;
        }
        block_size = proc_rows / mpi_mcw_size;
        for (i = 0; i < mpi_mcw_size - 1; ++i) {
            l = curr_row + 1 + block_size * i;
            MPI_Bcast((void *) (matrix + l * COLUMNS), block_size * COLUMNS,
                      MPI_DOUBLE, i, COMM);
        }
        
        l = curr_row + 1 + block_size * i;
        block_size += proc_rows % mpi_mcw_size;
        MPI_Bcast((void *) (matrix + l * COLUMNS), block_size * COLUMNS,
                  MPI_DOUBLE, i, COMM);
    }
}


/**
 * \brief Finding values of unknown variables.
 *
 * \param matrix is pointer to the matrix.
 * \param ROWS is number of rows in the matrix part.
 * \param COLUMNS is number of columns in the matrix.
 * 
 * \return vector of found values.
 */
double *MPI_Gauss_Backward(double *matrix, int ROWS, int COLUMNS)
{
    int i, j;
    double *results;
    
    results = (double *) calloc((size_t) ROWS, ROWS * sizeof(double));
    
    for (i = ROWS - 1; i >= 0; --i) {
        if (mpi_mcw_rank == 0) {
            if (Failure_gen(CHANCE)) {
                printf("(Backward) N %d \t STEP %d\n",
                       mpi_mcw_size - 1, i);
                raise(SIGKILL);
            }
        }
        
        for (j = 1; j < ROWS - i; ++j) {
            results[i] -= matrix[(i + 1) * COLUMNS - j - 1] *
                results[COLUMNS - j - 1];
        }
        results[i] += matrix[(i + 1) * COLUMNS - 1];
    }
    
    return results;
}


/**
 * \brief Gaussian elimination.
 *
 * Format input file:
 * <Rows> <Columns>
 * <Matrix content>
 *
 * Format output file:
 * <Number>
 * <Vector of found values>
 *
 * \param input is input filename.
 * \param output is output filename.
 *
 * \return zero if success.
 */
int MPI_Gauss(const char *input, const char *output)
{
    FILE *fp;
    double *matrix = NULL, *results = NULL;
    int columns, rows, i;

    COMM = MPI_COMM_WORLD;

    if (mpi_mcw_rank == 0) {
        fp = fopen(input, "r");
    }
    MPI_Matrix_Partition(fp, &matrix, &rows, &columns);
    if (mpi_mcw_rank == 0) {
        fclose(fp);
    }
    
    MPI_Gauss_Forward(matrix, rows, columns);

    results = MPI_Gauss_Backward(matrix, rows, columns);

    MPIX_Comm_shrink(COMM, &NEWCOMM);
    COMM = NEWCOMM;
    MPI_Comm_rank(COMM, &mpi_mcw_rank);
    MPI_Comm_size(COMM, &mpi_mcw_size);
    
    printf("SIZE %d\n", mpi_mcw_size);
    if (mpi_mcw_rank == 0) {
        fp = fopen(output, "w");
        fprintf(fp, "%d\n", rows);
        for (i = 0; i < rows; ++i) {
            fprintf(fp, "%.6f ", results[i]);
        }
        fclose(fp);
    }
    
    free(results);
    free(matrix);
    
    return 0;
}



/**
 * \brief Main function. Calculate matrix from file.
 *
 * \param arg_1 is input filename.
 * \param arg_2 is output filename.
 *
 * \return zero if success.
 */
int main(int argc, char *argv[])
{    
    if (argc < 3) {
        printf("Too few arguments!\n");
        return 1;
    }
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_mcw_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_mcw_rank);
    MPI_Comm_create_errhandler(mpi_error_handler, &new_eh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, new_eh);
    
    /* Root process */
    if (mpi_mcw_rank == 0) {
        /* Start timer */
        gettimeofday(&start, NULL);
    }

    /* MPI_Gauss(argv[1], argv[2]); */
    MPI_Gauss(argv[1], argv[2]);

    /* Root process */
    if (mpi_mcw_rank == 0) {
        /* Stop timer */
        gettimeofday(&end, NULL);
        
        /* Human readable format */
        end.tv_sec -= start.tv_sec;
        end.tv_usec += end.tv_sec * 1000000;
        end.tv_usec -= start.tv_usec;
        end.tv_sec = end.tv_usec / 1000000;
        end.tv_usec -= end.tv_sec * 1000000;
        printf("time: %ld.%ld sec.\n", end.tv_sec, end.tv_usec);
    }
    
    MPI_Finalize();

    return 0;
}
