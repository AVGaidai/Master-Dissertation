/** \file */
/*
 * GAUSS ELIMINATION (MPI VERSION) BY ANATOLY GAIDAI
 * 2018
 */
#include <stdio.h>
#include <stdlib.h>   /* malloc() */

#include <unistd.h>   /* sleep()  */

#include <string.h>   /* memcpy() */

#include <sys/time.h> /* gettimeofday () */

#include <mpi.h>


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
    int rank, commsize, i, j;
    double *ptr;
    
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Root process */
    if (rank == 0) {
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

    MPI_Bcast((void *) COLUMNS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *) ROWS, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        *matrix = (double *) malloc((*COLUMNS) * (*ROWS) *sizeof(double));
    }
    
    MPI_Bcast((void *) (*matrix), (*COLUMNS) * (*ROWS),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
    int rank, commsize;
    double k;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    for ( ; curr_row < ROWS; ++curr_row) {
        k = matrix[curr_row * COLUMNS + curr_row];
        for (i = curr_row + 1; i < COLUMNS; ++i) {
            matrix[curr_row * COLUMNS + i] /= k;
        }
        matrix[curr_row * COLUMNS + curr_row] = 1.00;

        if (curr_row == ROWS - 1) continue;
        
        proc_rows = ROWS - (curr_row + 1);
        block_size = proc_rows / commsize;
        l = curr_row + 1 + block_size * rank;
        if (rank == commsize - 1) {
            block_size += proc_rows % commsize;
        }
    
        for (i = 0; i < block_size; ++i) {
            k = 0.00 - matrix[(l + i) * COLUMNS + curr_row]; 
            for (j = curr_row + 1; j < COLUMNS; ++j) {
                matrix[(l + i) * COLUMNS + j] +=
                    matrix[curr_row * COLUMNS + j] * k;
            }
            matrix[(l + i) * COLUMNS + curr_row] = 0.00;
        }
        block_size = proc_rows / commsize;
        for (i = 0; i < commsize - 1; ++i) {
            l = curr_row + 1 + block_size * i;
            MPI_Bcast((void *) (matrix + l * COLUMNS), block_size * COLUMNS,
                      MPI_DOUBLE, i, MPI_COMM_WORLD);
        }
        l = curr_row + 1 + block_size * i;
        block_size += proc_rows % commsize;
        MPI_Bcast((void *) (matrix + l * COLUMNS), block_size * COLUMNS,
                  MPI_DOUBLE, i, MPI_COMM_WORLD);
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
    int rank, commsize;
    FILE *fp;
    double *matrix = NULL, *results = NULL;
    int columns, rows, i;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    if (rank == 0) {
        fp = fopen(input, "r");
    }
    MPI_Matrix_Partition(fp, &matrix, &rows, &columns);
    if (rank == 0) {
        fclose(fp);
    }
    
    MPI_Gauss_Forward(matrix, rows, columns);

    results = MPI_Gauss_Backward(matrix, rows, columns);

    if (rank == 0) {
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
    int rank;
    struct timeval start, end;
    
    if (argc < 3) {
        printf("Too few arguments!\n");
        return 1;
    }
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Root process */
    if (rank == 0) {
        /* Start timer */
        gettimeofday(&start, NULL);
    }
    
    /* MPI_Gauss(argv[1], argv[2]); */
    MPI_Gauss(argv[1], argv[2]);

    /* Root process */
    if (rank == 0) {
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
