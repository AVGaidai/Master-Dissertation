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
 * \param file descriptor of the input matrix.
 * \param pointer to the output part.
 * \param pointer to the number of rows in the output part.
 * \param pointer to the number of columns in the output part.
 */
void MPI_Matrix_Partition(FILE *fp, double **part, int *ROWS, int *COLUMNS)
{
   /*
    * cnt is counter of rows in the output part
    * row is currently processed row
    */
    int rank, commsize, i, j, cnt, receiver;
    double *row;

    MPI_Status status;
        
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello World, I'am %d process!\n", rank);
    *part = NULL;
    cnt = 0;
    /* Main process */
    if (rank == 0) {
        /* Reading matrix size */
        fscanf(fp, "%d %d", ROWS, COLUMNS);
        /* Allocate memory for row */
        row = (double *) malloc(*COLUMNS * sizeof(double));
        /* Alternate reading matrix rows */
        for (i = 0; i < *ROWS; ++i) {
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
                memcpy((void *) *part + cnt * (*COLUMNS),
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

    /* Broadcast message with nomber of columns */
    MPI_Bcast((void *) COLUMNS, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
            memcpy((void *) *part + cnt * (*COLUMNS),
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
 * \param matrix is pointer on matrix.
 * \param SIZE_X is number of rows in the matrix.
 * \param SIZE_Y is number of columns in the matrix.
 */
void MPI_Gauss_Forward(double *matrix, int SIZE_X, int SIZE_Y)
{
    int commsize, rank, current_row = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello! I'm %d from %d.\n", rank, commsize);
}


/**
 * \brief Finding values of unknown variables.
 *
 * \param matrix is pointer on matrix.
 * \param SIZE_X is number of rows in the matrix.
 * \param SIZE_Y is number of columns in the matrix.
 *
 * \return Result is vector values.
 */
double *calculate_matrix(double *matrix, int SIZE_X, int SIZE_Y)
{
    ;
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
    int commsize, rank, ROWS, COLUMNS;
    double *matrix = NULL; 
    FILE *fp;
    
    if (argc < 2) return 1;
            
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
        fp = fopen(argv[1], "r");
    
    MPI_Matrix_Partition(fp, &matrix, &ROWS, &COLUMNS);

    if (rank == 0)
        fclose(fp);
    
    sleep(rank);
    print_matrix(matrix, ROWS, COLUMNS);

    if (ROWS) free(matrix);

    MPI_Finalize();
}
