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


void MPI_Matrix_Partition(FILE *fp, double **part, int *ROWS, int *COLUMNS)
{
    int rank, commsize, i, j, cnt, receiver;
    double *row;

    MPI_Status status;
        
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello World, I'am %d process!\n", rank);
    *part = NULL;
    cnt = 0;
    if (rank == 0) {
        fscanf(fp, "%d %d", ROWS, COLUMNS);
        row = (double *) malloc(*COLUMNS * sizeof(double));
        for (i = 0; i < *ROWS; ++i) {
            for (j = 0; j < *COLUMNS; ++j) {
                fscanf(fp, "%lf", &row[j]);
            }
            receiver = i % commsize;
            if (receiver == 0) {
                *part = (double *) realloc(*part, (cnt + 1) * sizeof(double));
                memcpy((void *) *part + cnt * (*COLUMNS),
                       (void *) row,
                       sizeof(double) * (*COLUMNS));  
                ++cnt;
            } else {
                //printf("send to %d:\n", receiver);
                //print_matrix(row, 1, *COLUMNS);
                MPI_Send((const void *) row,
                         *COLUMNS,
                         MPI_DOUBLE,
                         receiver,
                         receiver,
                         MPI_COMM_WORLD);
            }
        }
        for (i = 1; i < commsize; ++i) { 
            MPI_Send((const void *) row,
                     *COLUMNS,
                     MPI_DOUBLE,
                     i,
                     0,
                     MPI_COMM_WORLD);
        }
    }
    
    MPI_Bcast((void *) COLUMNS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank != 0) {
        //sleep(2);
        row = (double *) malloc(*COLUMNS * sizeof(double));
        
        while (1) {
            MPI_Recv((void *) row,
                     *COLUMNS,
                     MPI_DOUBLE,
                     0,
                     MPI_ANY_TAG,
                     MPI_COMM_WORLD,
                     &status);

             if (!status.MPI_TAG) break;
            //sleep(rank);
            //printf("%d\n", rank);
            //print_matrix(row, 1, *COLUMNS);
            *part = (double *) realloc(*part, (cnt + 1) * sizeof(double));
            memcpy((void *) *part + cnt * (*COLUMNS),
                   (void *) row,
                   sizeof(double) * (*COLUMNS));  
            ++cnt;
        }
    }
    *ROWS = cnt;
    //sleep(rank);
    //print_matrix(*part, *ROWS, *COLUMNS);
    MPI_Barrier(MPI_COMM_WORLD);    
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
    //MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        fclose(fp);
    
    sleep(rank);
    print_matrix(matrix, ROWS, COLUMNS);

    //if (ROWS) free(matrix);
    MPI_Finalize();
}
