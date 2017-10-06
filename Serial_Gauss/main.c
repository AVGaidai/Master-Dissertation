/** \file */
/*
 * GAUSS ELIMINATION (SERIAL VERSION) BY ANATOLY GAIDAI
 * 2017
 */
#include <stdio.h>
#include <stdlib.h>


/**
 * \brief Print content of matrix 'matrix'.
 * 
 * \param matrix is pointer on matrix.
 * \param ROWS is number of rows in the matrix.
 * \param COLUMNS is number of columns in the matrix.
 *
 * \return no return value.
 */
void print_matrix(double *matrix, int ROWS, int COLUMNS)
{
   /*
    * i is current row in the matrix
    * j is current columns in the matrix
    * offset is shift in the current row of matrix
    */
    int i, j, offset;
    
    printf("Matrix:\n");
    for (i = 0; i < ROWS; ++i) {
        offset = i * COLUMNS;
	for (j = 0; j < COLUMNS; ++j) {
	    printf("%.2f\t", matrix[offset + j]);
	}
	printf("\n");
    }
}


/**
 * \brief Matrix transformation to lower-triangular.
 *
 * \param matrix is pointer on matrix.
 * \param ROWS is number of rows in the matrix.
 * \param COLUMNS is number of columns in the matrix.
 *
 * \return no return value.
 */
void Gauss_Forward(double *matrix, int ROWS, int COLUMNS)
{    
   /*
    * i is current row in the matrix
    * j is current column in the matrix
    * k is current row in the calculating part of matrix
    * l is current column in the calculating part of matrix
    * offset is shift in the current row of matrix
    * coeff_1 needed to create a single diagonal
    * coeff_2 needed to transform to lower-triangular
    */
    int i, j, k, l, offset;
    double coeff_1, coeff_2;

    /* Alternate processing (forward) of the matrix */
    for (i = 0; i < ROWS; ++i) {
        offset = i * COLUMNS;
        /* Finding coefficient for current row */
	coeff_1 = 1.0 / matrix[offset + i];
        /* Transforming current row */
	for (j = i; j < COLUMNS; ++j) {
	    matrix[offset + j] *= coeff_1;
	}
        /* Transforming other rows */
	for (k = i + 1; k < ROWS; ++k) {
            offset = k * COLUMNS;
            /* Finding coefficient for each row */
	    coeff_2 = -1.0 * matrix[offset + i];
            /* Transforming row */
	    for (l = i; l < COLUMNS; ++l) {
		matrix[offset + l] += matrix[i * COLUMNS + l] * coeff_2;
	    }
	}
    }
}


/**
 * \brief Finding values of unknown variables.
 *
 * \param matrix is pointer on matrix.
 * \param ROWS is number of rows in the matrix.
 * \param COLUMNS is number of columns in the matrix.
 *
 * \return vector values.
 */
double *Gauss_Backward(double *matrix, int ROWS, int COLUMNS)
{
   /*
    * i is current row in the matrix
    * j is current column in the matrix
    * offset is shift in the current row of matrix
    * result is vector of found values
    */    
    int i, j, offset;
    double *result;

    /* Allocate memory for result vector */
    result = (double *) malloc(ROWS * sizeof (double));

   /*
    * Vector is formed from right to left 
    */
    for (i = ROWS - 1; i >= 0; --i) {
        offset = (i + 1) * COLUMNS;
	result[COLUMNS - (ROWS - i) - 1] = matrix[offset - 1];
        for (j = 1; j < ROWS - i; ++j) {
	    result[COLUMNS - (ROWS - i) - 1] -= result[COLUMNS - j - 1] *
                                                matrix[offset - 1 - j];
	}
    }	    
    
    return result; 
}

/**
 * \brief Gauss elimination.
 *
 * Format input file:
 * <Rows> <Columns>
 * <Matrix content>
 *
 * Format output file:
 * <Number>
 * <Vector of found values>
 *
 * \param input is name of input file.
 * \param output is name of output file.
 *
 * \return zero if success.
 */
int Gauss(const char *input, const char *output)
{
   /*
    * fp is file descriptor for I/O
    * ROWS is number of rows in the matrix
    * COLUMNS is number of columns in the matrix
    * i is current row number in the matrix
    * j is current column number in the matrix
    * matrix is data matrix
    * result is vector of found values
    */
    FILE *fp;
    int ROWS, COLUMNS, i, j;
    double *matrix, *result;

    /* Open file for reading data */
    fp = fopen(input, "rb");
    
    fscanf(fp, "%d", &ROWS);
    fscanf(fp, "%d", &COLUMNS);

    /* Allocate memory for matrix */
    matrix = (double *) malloc(ROWS * COLUMNS * sizeof (double));
    /* Reading matrix from the input file */ 
    for (i = 0; i < ROWS; ++i) {
	for (j = 0; j < COLUMNS; ++j) {
	    fscanf(fp, "%lf", &matrix[i * COLUMNS + j]);
	}
    }
    
    fclose(fp); /* Close file descriptor */
    
    /* Forwadr elimination */
    Gauss_Forward(matrix, ROWS, COLUMNS);
   
    /* Backward elimenation */
    result = Gauss_Backward(matrix, ROWS, COLUMNS);

    /* Open file for writing found values */
    fp = fopen(output, "wb");
    /* Writing number of found values */
    fprintf(fp, "%d\n", ROWS);
    /* Writing vector of found values */
    for (i = 0; i < ROWS; ++i) {
        fprintf(fp, "%f ", result[i]);
    }
    
    free(result);
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
    /* Too few arguments */
    if (argc < 3) {
        printf("Too few argements!\n");
        return 1;
    }
    
    if (Gauss(argv[1], argv[2]))
        printf("failed\n");

    return 0;
}
