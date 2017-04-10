/** \file */
#include <stdio.h>
#include <stdlib.h>


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
   /*
    * coeff_1 needed to create a single diagonal
    * coeff_2 needed to transform to lower-triangular
    */
    double coeff_1, coeff_2;
    
    for (int i = 0; i < SIZE_X; ++i) {
	//print_matrix(matrix, SIZE_X, SIZE_Y);
	coeff_1 = 1.0 / matrix[i * SIZE_Y + i];
	for (int j = i; j < SIZE_Y; ++j) {
	    matrix[i * SIZE_Y + j] *= coeff_1;
	}
	
	MPI_Scatter((const void *) matrix + i * SIZE_Y, SIZE_Y, MPI_DOUBLE,
		    (void*) recv_buf, SIZE_Y, 
	for (int k = i + 1; k < SIZE_X; ++k) {
	    coeff_2 = -1.0 * matrix[k * SIZE_Y + i];
	    for (int l = i; l < SIZE_Y; ++l) {
		matrix[k * SIZE_Y + l] += matrix[i * SIZE_Y + l] * coeff_2;
	    }
	}
    }
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

    int comsize, rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_transform_matrix(matrix, SIZE_X, SIZE_Y);

    MPI_Finalize();

    print_matrix(matrix, SIZE_X, SIZE_Y);
    
    double *result;
    
    result = calculate_matrix(matrix, SIZE_X, SIZE_Y);
    print_vector(result, SIZE_Y - 1);
    
    free(matrix);
}
