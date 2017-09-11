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
 * \brief Matrix transformation to lower-triangular.
 *
 * \param matrix is pointer on matrix.
 * \param SIZE_X is number of rows in the matrix.
 * \param SIZE_Y is number of columns in the matrix.
 */
void MPI_transform_matrix(double *matrix, int SIZE_X, int SIZE_Y)
{
    int commsize, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello! I'm %d from %d.\n", rank, commsize);
   /*
    * coeff_1 needed to create a single diagonal
    * coeff_2 needed to transform to lower-triangular
    */
    double coeff_1, coeff_2;
    
    double *buf_1, *buf_2, *buf_3;
    int j, k, l, block_size, remainder, all;
    MPI_Status status;

    buf_1 = (double *) malloc(sizeof(double) * SIZE_Y);
    for (int i = 0; i < SIZE_X; ++i) {
	block_size = (SIZE_X - (i + 1)) / commsize;
	if (block_size == 0) {
	    remainder = SIZE_X - (i + 1);
	} else {
	    remainder = (SIZE_X - (i + 1)) % commsize + block_size;
	}
	if (rank == 0) {
	    coeff_1 = 1.0 / matrix[i * SIZE_Y + i];
	    for (int j = i; j < SIZE_Y; ++j) {
		matrix[i * SIZE_Y + j] *= coeff_1;
	    }
	    if (i == SIZE_X - 1) break;
	    memcpy(buf_1, matrix + i * SIZE_Y, SIZE_Y * sizeof(double));

     	    k = 0;
	    all = 0;
	    buf_3 = (double *) malloc(sizeof(double) * SIZE_Y * block_size);
	    for (j = 0; j < commsize - 1; ++j) {
		if (i + 1 + j * block_size >= SIZE_X - remainder) break;
		memcpy(buf_3, matrix + (i + 1) * SIZE_Y + j * block_size *
                       SIZE_Y, sizeof(double) * SIZE_Y * block_size);
		MPI_Send(buf_3, block_size * SIZE_Y,
                         MPI_DOUBLE, j + 1, all++, MPI_COMM_WORLD);
		printf("to %d:\n", j + 1);
		print_matrix(matrix + (i + 1) * SIZE_Y + j * block_size * SIZE_Y,
                             block_size, SIZE_Y);
	    }
	    free(buf_3);
	    for (j = 1; j < commsize; ++j) {
		MPI_Send((const void *) matrix, block_size * SIZE_Y,
			 MPI_DOUBLE, j, commsize + 1, MPI_COMM_WORLD);
	    }
	}
	if (i == SIZE_X - 1) break;
	
	MPI_Bcast(buf_1, SIZE_Y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
	    printf("Bcast string:\n");
	    print_vector(buf_1, SIZE_Y);
	}
	
	if (rank != 0) {
	    buf_2 = (double *) malloc(sizeof(double) * block_size * SIZE_Y);
	    printf("Process %d\n", rank);
	    do {
		printf("buf_1:\n");
		print_vector(buf_1, SIZE_Y);
		MPI_Recv((void *) buf_2, block_size * SIZE_Y, MPI_DOUBLE,
			 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		printf("status %d\n", status.MPI_TAG);
		if (status.MPI_TAG == commsize + 1) break;
		printf("part start:\n");
		print_matrix(buf_2, block_size, SIZE_Y);
		for (k = 0; (buf_1[k] - 1.00 > 0.000001) ||
                         (buf_1[k] - 1.00 < -0.000001); ++k);
		printf("k=%d; dif=%lf\n", k, buf_1[k] - 1.00);
		for (l = 0; l < block_size; ++l) { 
		    coeff_2 = -1.0 * buf_2[l * SIZE_Y + k];
		    for (j = k; j < SIZE_Y; ++j) {
			buf_2[l * SIZE_Y + j] += buf_1[j] * coeff_2;
		    }
		}
		printf("part final:\n");
		print_matrix(buf_2, block_size, SIZE_Y);
		MPI_Send((const void *) buf_2, block_size * SIZE_Y, MPI_DOUBLE,
			 0, status.MPI_TAG, MPI_COMM_WORLD);
	    } while (status.MPI_TAG != commsize + 1);
	    free(buf_2);
	} else {
	    printf("Process %d\n", rank);
	    printf("buf_1:\n");
	    print_vector(buf_1, SIZE_Y);	    
	    buf_2 = (double *) malloc(sizeof(double) * remainder * SIZE_Y);
	    
	    memcpy(buf_2, matrix + (SIZE_X - remainder) * SIZE_Y,
		   remainder * SIZE_Y * sizeof(double));
	    printf("part start:\n");
	    print_matrix(buf_2, remainder, SIZE_Y);
	    for (k = 0; (buf_1[k] - 1.00 > 0.000001) ||
                     (buf_1[k] - 1.00 < -0.000001); ++k);
	    printf("k=%d; dif=%lf\n", k, buf_1[k] - 1.00);
	    for (l = 0; l < remainder; ++l) { 
		coeff_2 = -1.0 * buf_2[l * SIZE_Y + k];
		for (j = k; j < SIZE_Y; ++j) {
		    buf_2[l * SIZE_Y + j] += buf_1[j] * coeff_2;
		}
	    }
	    memcpy(matrix + (SIZE_X - remainder) * SIZE_Y, buf_2,
		   remainder * SIZE_Y * sizeof(double));
	    printf("part final:\n");
	    print_matrix(buf_2, remainder, SIZE_Y);
	    free(buf_2);

	    buf_3 = (double *) malloc(sizeof(double) * block_size * SIZE_Y);
	    printf("all=%d\n", all);
	    for (j = 0; j < all; ++j) {
		MPI_Recv((void *) buf_3, block_size * SIZE_Y, MPI_DOUBLE,
                         MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (status.MPI_TAG == commsize + 1) status.MPI_TAG = 0;
		memcpy(matrix + (i + 1) * SIZE_Y + status.MPI_TAG * block_size *
                       SIZE_Y, buf_3, block_size * SIZE_Y * sizeof(double));
	    }
	    free(buf_3);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
	    printf("=========================\n");
	    print_matrix(matrix, SIZE_X, SIZE_Y);
	}
    }
    free(buf_1);
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

    int commsize, rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *matrix = NULL;
    int SIZE_X, SIZE_Y;
    
    if (rank == 0) {
	FILE *fp;
	
	fp = fopen(argv[1], "rb");

	fscanf(fp, "%d", &SIZE_X);
	fscanf(fp, "%d", &SIZE_Y);
    
	matrix = (double *) malloc(SIZE_X * SIZE_Y * sizeof (double));
	for (int i = 0; i < SIZE_X; ++i) {
	    for (int j = 0; j < SIZE_Y; ++j) {
		fscanf(fp, "%lf", &matrix[i * SIZE_Y + j]);
	    }
	}
	fclose(fp);
	print_matrix(matrix, SIZE_X, SIZE_Y);
    }

    MPI_Bcast(&SIZE_X, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&SIZE_Y, 1, MPI_INT, 0, MPI_COMM_WORLD);
   
    MPI_transform_matrix(matrix, SIZE_X, SIZE_Y);

    if (rank == 0) {
	printf("4\n");    
	print_matrix(matrix, SIZE_X, SIZE_Y);
	printf("4\n");    
	double *result;
    
	result = calculate_matrix(matrix, SIZE_X, SIZE_Y);
	print_vector(result, SIZE_Y - 1);
	free(matrix);
    }

    MPI_Finalize();
}
