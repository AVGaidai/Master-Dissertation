/** \file */
#include <stdio.h>
#include <stdlib.h>

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
 * \brief Matrix transformation to lower-triangular.
 *
 * \param matrix is pointer on matrix.
 * \param SIZE_X is number of rows in the matrix.
 * \param SIZE_Y is number of columns in the matrix.
 */
void MPI_transform_matrix(double *matrix, int SIZE_X, int SIZE_Y)
{
    int commsize, rank;
    int argc;
    char **argv;

    printf("123\n");
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello! I'm %d from %d.\n", rank, commsize);
   /*
    * coeff_1 needed to create a single diagonal
    * coeff_2 needed to transform to lower-triangular
    */
    double coeff_1, coeff_2;
    
    double *buf_1, *buf_2, *buf_3;
    int j, k, l;
    MPI_Status status;
    
    for (int i = 0; i < SIZE_X; ++i) {
	if (rank == 0) {
	    coeff_1 = 1.0 / matrix[i * SIZE_Y + i];
	    for (int j = i; j < SIZE_Y; ++j) {
		matrix[i * SIZE_Y + j] *= coeff_1;
	    }
	    MPI_Bcast((void *) matrix + i * SIZE_Y, SIZE_Y, MPI_DOUBLE,
		      0, MPI_COMM_WORLD);
	    k = 0;
	    l = 0;
	    for (j = i + 1; j < SIZE_X; ++j) {
		MPI_Send((const void *) matrix + j * SIZE_Y, SIZE_Y, MPI_DOUBLE,
			 (k++ % commsize) + 1, 0, MPI_COMM_WORLD);
		++l;
	    }

	    for (j = 0; j < commsize; ++j) {
		MPI_Send((const void *) matrix, SIZE_Y,
			 MPI_DOUBLE, k, 1, MPI_COMM_WORLD);
		++l;
	    }
	    
	    buf_3 = (double *) malloc(sizeof(double) * SIZE_Y);
	    for (j = 0; j < l; ++j) {
		MPI_Recv((void *) buf_3, SIZE_Y, MPI_DOUBLE, MPI_ANY_SOURCE,
			 0, MPI_COMM_WORLD, &status);
		if (status.MPI_TAG != 0) {
		    memcpy(matrix + status.MPI_TAG * SIZE_Y,
			   buf_3, SIZE_Y * sizeof(double));
		}
	    }
	}
		    
	if (rank != 0) {
	    buf_1 = (double *) malloc(sizeof(double) * SIZE_Y);
	    buf_2 = (double *) malloc(sizeof(double) * SIZE_Y);
	    MPI_Recv((void *) buf_1, SIZE_Y, MPI_DOUBLE, 0, 0,
		     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    do {
		MPI_Recv((void *) buf_2, SIZE_Y, MPI_DOUBLE,
			 0, 0, MPI_COMM_WORLD, &status);
		for (k = 0; buf_1[k] - 1 < 0.00001; ++k);
		coeff_2 = -1.0 * buf_2[k];
		for (j = k; j < SIZE_Y; ++j) {
		    buf_2[j] += buf_1[j] * coeff_2;
		}
		MPI_Send((const void *) buf_2, SIZE_Y, MPI_DOUBLE,
			 0, status.MPI_TAG, MPI_COMM_WORLD);
	    } while (status.MPI_TAG != 0);
	}
    }

//    MPI_Finalize();
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
    double *result;

    result = (double *) malloc((SIZE_Y - 1) * sizeof (double));

   /*
    * Vector is formed from right to left 
    */
    for (int i = SIZE_X - 1; i >= 0; --i) {
	result[SIZE_Y - (SIZE_X - i) - 1] = matrix[(i + 1) * SIZE_Y - 1];
	int j = 1;
	
	for (j = 1; j < SIZE_X - i; ++j) {
	    result[SIZE_Y - (SIZE_X - i) - 1] -= result[SIZE_Y - j - 1] *
		                               matrix[(i + 1) * SIZE_Y - 1 - j];
	}
    }	    
    
    return result; 
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
    if (argc < 2) return 1;
    
    FILE *fp;

    fp = fopen(argv[1], "rb");

    int SIZE_X, SIZE_Y;
    
    fscanf(fp, "%d", &SIZE_X);
    fscanf(fp, "%d", &SIZE_Y);
    
    double *matrix;

    matrix = (double *) malloc(SIZE_X * SIZE_Y * sizeof (double));
    for (int i = 0; i < SIZE_X; ++i) {
	for (int j = 0; j < SIZE_Y; ++j) {
	    fscanf(fp, "%lf", &matrix[i * SIZE_Y + j]);
	}
    }
    fclose(fp);
    print_matrix(matrix, SIZE_X, SIZE_Y);

//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("4\n");
    MPI_transform_matrix(matrix, SIZE_X, SIZE_Y);

    print_matrix(matrix, SIZE_X, SIZE_Y);
    
    double *result;
    
    result = calculate_matrix(matrix, SIZE_X, SIZE_Y);
    print_vector(result, SIZE_Y - 1);
    
    free(matrix);

    MPI_Finalize();
}
