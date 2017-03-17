/** \file */
#include <stdio.h>
#include <stdlib.h>

/*
typedef struct {

    int dividend;
    int divisor;

    double division;
    
} Quotient;


Quotient quot_init(int dividend, int divisor)
{
    return (Quotient) {
	.dividend = dividend,
	.divisor = divisor,
	.division = dividend / (double) divisor
	}
}


Quotient quot_add(Quotient A, Quotient B)
{
    return (Quotient) {
	.dividend = A.dividend + B.dividend,
	.divisor = A.divisor + B.divisor,
	.division = .dividend / (double) .divisor
	}
}


Quotient quot_mul(Quotient A, Quotient B)
{
    return (Quotient) {
	.dividend = A.dividend * B.dividend,
	.divisor = A.divisor * B.divisor,
	.divisiob = .dividend / (double) divisor
	}
}
*/


/**
 * Print content of vecrot 'vec'.
 * \param *vec -- pointer on vecor;
 * \param size -- size of vector.
 *
 * \return -- no return value.
 */
void print_vector(double *vec, int size)
{
    printf("Vector:\n");
    for (int i = 0; i < size; ++i)
	printf("x%d=%.4lf\n", i, vec[i]);
}


/**
 * Print content of matrix 'matrix'.
 * \param *matrix -- pointer on matrix;
 * \param SIZE_X -- number of rows in the matrix;
 * \param SIZE_Y -- number of columns in the matrix.
 *
 * \return -- no return value.
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
 * Matrix transformation to lower-triangular.
 * \param *matrix -- pointer on matrix;
 * \param SIZE_X -- number of rows in the matrix;
 * \param SIZE_Y -- number of columns in the matrix.
 *
 * \return -- no return value.
 */
void transform_matrix(double *matrix, int SIZE_X, int SIZE_Y)
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
	for (int k = i + 1; k < SIZE_X; ++k) {
	    coeff_2 = -1.0 * matrix[k * SIZE_Y + i];
	    for (int l = i; l < SIZE_Y; ++l) {
		matrix[k * SIZE_Y + l] += matrix[i * SIZE_Y + l] * coeff_2;
	    }
	}
    }
}


/**
 * Finding values of unknown variables.
 * \param *matrix -- pointer on matrix;
 * \param SIZE_X -- number of rows in the matrix;
 * \param SIZE_Y -- number of columns in the matrix.
 *
 * \return -- vector values.
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
 * Main function. Calculate matrix from file.
 * Format file:
 * <Rows> <Columns>
 * <Matrix content>
 * \param arg_1 -- filename.
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
    transform_matrix(matrix, SIZE_X, SIZE_Y);
    print_matrix(matrix, SIZE_X, SIZE_Y);
    
    double *result;
    
    result = calculate_matrix(matrix, SIZE_X, SIZE_Y);
    print_vector(result, SIZE_Y - 1);
    
    free(matrix);
}
