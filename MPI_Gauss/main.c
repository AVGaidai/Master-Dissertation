/** \file */
#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>   // sleep

#include <string.h>

#include <mpi.h>

/**
 * \brief Print content of vecrot 'vec'.
 *
 * \param vec is pointer on vecor.
 * \param size is size of vector.
 */
void print_vector(double *vec, int size)
{
    printf("Vector:\n");
    for (int i = 0; i < size; ++i)
	printf("x%d=%.4lf\n", i, vec[i]);
}


/**
 * \brief Print content of matrix 'matrix'.
 *
 * \param matrix is pointer on matrix.
 * \param SIZE_X is number of rows in the matrix.
 * \param SIZE_Y is number of columns in the matrix.
 */
void print_matrix(double *matrix, int SIZE_X, int SIZE_Y)
{
    printf("Matrix:\n");
    for (int i = 0; i < SIZE_X; ++i) {
	for (int j = 0; j < SIZE_Y; ++j) {
	    printf("%.2lf\t", matrix[i * SIZE_Y + j]);
	}
	printf("\n");
    }
}

/**
 * \brief Partition of the input matrix into parts.
 *
 * \param fp is file descriptor of the input matrix.
 * \param part is pointer to the output part.
 * \param ROWS is pointer to the number of rows in the output part.
 * \param COLUMNS is pointer to the number of columns in the output part.
 * \param ALLROWS is pointer to the number of rows in the input matrix.
 */
void MPI_Matrix_Partition(FILE *fp, double **part,
                          int *ROWS, int *COLUMNS, int *ALLROWS)
{
   /*
    * i is index of the row in the input matrix
    * j is index of the columns in the row
    * cnt is counter of rows in the output part
    * row is currently processed row
    */
    int rank, commsize, i, j, cnt, receiver;
    double *row;

    MPI_Status status;
        
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    *part = NULL;
    cnt = 0;
    /* Main process */
    if (rank == 0) {
        /* Reading matrix size */
        fscanf(fp, "%d %d", ALLROWS, COLUMNS);
        /* Allocate memory for row */
        row = (double *) malloc(*COLUMNS * sizeof(double));
        /* Alternate reading matrix rows */
        for (i = 0; i < *ALLROWS; ++i) {
            for (j = 0; j < *COLUMNS; ++j) {
                fscanf(fp, "%lf", &row[j]);
            }
            /* Define receiver */
            receiver = i % commsize;
            /* For main process */
            if (receiver == 0) {
                /* Reallocate memory and add new row in the output part */
                *part = (double *) realloc(*part, (cnt + 1) *
                                           (*COLUMNS) * sizeof(double));
                memcpy((void *) (*part + cnt * (*COLUMNS)),
                       (void *) row,
                       sizeof(double) * (*COLUMNS));  
                ++cnt; // Increase counter of the rows
            /* For other rows */
            } else {
                /* Send the row to receiver */
                MPI_Send((const void *) row, // Row
                         *COLUMNS,           // Number of columns in the row
                         MPI_DOUBLE,         // Elements type
                         receiver,           // Receiver
                         receiver,           // Tag
                         MPI_COMM_WORLD);    // Commutator
            }
        }
        /* Sending message about the end reading to all process */ 
        for (i = 1; i < commsize; ++i) { 
            MPI_Send((const void *) row, // Any row (dummy)
                     *COLUMNS,           // Number of columns in the row
                     MPI_DOUBLE,         // Elements type
                     i,                  // Receiver
                     0,                  // Tag
                     MPI_COMM_WORLD);    // Commutator
        }
    }

    /* Broadcast message with number of columns */
    MPI_Bcast((void *) COLUMNS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* Broadcast message with number of rows in the input matrix */    
    MPI_Bcast((void *) ALLROWS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    /* For other process */
    if (rank != 0) {
        /* Allocate memory for row */
        row = (double *) malloc(*COLUMNS * sizeof(double));
        /* Waiting for message */
        while (1) {
            /* Receive row of the matrix */
            MPI_Recv((void *) row,   // Row
                     *COLUMNS,       // Number of columns in the row
                     MPI_DOUBLE,     // Elements type
                     0,              // Sender
                     MPI_ANY_TAG,    // Tag
                     MPI_COMM_WORLD, // Commetator
                     &status);       // Status

            /* If receive message about the end reading */ 
            if (!status.MPI_TAG) break;
            /* Reallocate memory and add new row in the output part */
            *part = (double *) realloc(*part, (cnt + 1) *
                                       (*COLUMNS) * sizeof(double));
            memcpy((void *) (*part + cnt * (*COLUMNS)),
                   (void *) row,
                   sizeof(double) * (*COLUMNS));  
            ++cnt; // Increase counter of the rows
        }
    }
    *ROWS = cnt;
    free(row);
}


/**
 * \brief Matrix transformation to lower-triangular.
 *
 * \param part is pointer on matrix part.
 * \param ROWS is number of rows in the matrix part.
 * \param COLUMNS is number of columns in the matrix part.
 * \param ALLROWS is number of rows in the input matrix.
 */
void MPI_Gauss_Forward(double *part, int ROWS, int COLUMNS, int ALLROWS)
{
   /*
    * i is current row number in the general matrix
    * j is active node number
    * k is current row number in the matrix part (for send)
    * l is current column number in the current row
    * m is current row number in the matrix part (for modification)
    */
    int commsize, rank, offset;
    int i, j, k, l, m;
    double *row;
        
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Allocate memory for row */
    row = (double *) malloc(COLUMNS * sizeof(double));
    k = 0;
    /* Cyclic processing (forward) of the general matrix */
    for (i = 0; i < ALLROWS; ) {
        for (j = 0; j < commsize && i < ALLROWS; ++j) {
            /* If node is active */
            if (rank == j) {
                offset = k * COLUMNS;
                /* Modification of the current row */
                for (l = i + 1; l < COLUMNS; ++l) {
                    part[offset + l] /= part[offset + i];
                }
                part[offset + i] = 1;
                /* Copying modificated row to the buffer */
                memcpy(row, part + offset, COLUMNS * sizeof(double));
                ++k;
            }
            /* Broad cast message with current row */
            MPI_Bcast((void *) row, COLUMNS, MPI_DOUBLE, j, MPI_COMM_WORLD);

            /* Processing (forward) of the part matrix remainder */
            for (m = k; m < ROWS; ++m) {
                offset = m * COLUMNS;
                for (l = i + 1; l < COLUMNS; ++l) {
                    part[offset + l] -= part[offset + i] * row[l];
                }
                part[offset + i] = 0;
            }
            ++i;
        }
    }

    free(row);
}


/**
 * \brief Finding values of unknown variables.
 *
 * \param part is pointer on matrix part.
 * \param ROWS is number of rows in the matrix part.
 * \param COLUMNS is number of columns in the matrix part.
 * \param ALLROWS is number of rows in the input matrix.
 * 
 * \return vector of found values.
 */
double *MPI_Gauss_Backward(double *part, int ROWS, int COLUMNS, int ALLROWS)
{
   /*
    * i is current row number in the general matrix
    * j is active node number
    * k is current row number in the matrix part (for final value)
    * l is current columns in the row of calculation
    * m is current row number in the matrix part (for modification)
    * X is vector of found values
    * tmp is temporary variable
    */
    int commsize, rank, offset;
    int i, j, k, l, m;
    double *X;
    double tmp;
        
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* If node was activated */
    if (ROWS)
        /* Allocate memory for vector of found values */
        X = (double *) malloc(ROWS * sizeof(double));
    k = ROWS - 1;
    /* Cyclic processing (backward) of the general matrix */
    for (i = ALLROWS - 1; i >= 0; ) {
        for (j = i % commsize; j >= 0 && i >= 0; --j) {
            /* If node is active */
            if (rank == j) {
                offset = k * COLUMNS;
                /* Calculate value */
                X[k] = part[offset + COLUMNS - 1];
                for (l = COLUMNS - 2; l > i; --l) {
                    X[k] -= part[offset + l];
                }
                tmp = X[k];
                --k;
            }
            /* Broad cast message with found value */
            MPI_Bcast((void *) &tmp, 1, MPI_DOUBLE, j, MPI_COMM_WORLD);

            /* Calculating processing (backward) of the part matrix remainder */
            for (m = k; m >= 0; --m) {
                offset = m * COLUMNS;
                part[offset + i] *= tmp;
            }
            --i;
        }
    }
    
    return X;
}



/**
 * \brief Main function. Calculate matrix from file.
 *
 * Format file:
 * <Rows> <Columns>
 * <Matrix content>
 *
 * \param arg_1 is filename.
 */
int main(int argc, char *argv[])
{
    int commsize, rank, ROWS, COLUMNS, ALLROWS, i;
    double *X, *matrix = NULL; 
    FILE *fp;
    
    if (argc < 2) return 1;
            
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
        fp = fopen(argv[1], "r");
    
    MPI_Matrix_Partition(fp, &matrix, &ROWS, &COLUMNS, &ALLROWS);

    if (rank == 0)
        fclose(fp);

    /* sleep(rank); */
    /* if (rank == 0) print_matrix(matrix, ROWS, COLUMNS); */

    /* sleep(5); */
    
    MPI_Gauss_Forward(matrix, ROWS, COLUMNS, ALLROWS);

    /* sleep(rank); */
    /* print_matrix(matrix, ROWS, COLUMNS); */

    X = (double *) malloc(ROWS * sizeof(double));
    X = MPI_Gauss_Backward(matrix, ROWS, COLUMNS, ALLROWS);

    /* sleep(rank); */

    for (i = 0; i < ROWS; ++i) {
        printf("X%d = %lf\n", rank, X[i]);
        rank += commsize;
    }
    
    if (ROWS) {
        free(X);
        free(matrix);
    }
    
    MPI_Finalize();
}
