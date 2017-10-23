/** \file */
/* Programm for generation matrix with values [0 10] */
#include <stdio.h>
#include <stdlib.h>

#include <time.h>



/**
 * \brief Main function. It is generating matrix with values [0 10].
 *
 * \param arg_1 is number of rows into matrix.
 * \param arg_2 is number of columns into matrix.
 * \param arg_3 is filename containig result matrix (optional argument).
 *
 * \return 1 if too few arguments, 0 if success.
 */ 
int main(int argc, char *argv[])
{
   /*
    * fp is outfile descriptor
    * ROWS is number of rows into matrix
    * COLUMNS is number of columns into matrix
    * i and j are indexes
    */
    FILE *fp;
    int ROWS, COLUMNS, i, j;
    
    if (argc < 3) {
        fprintf(stderr, "Too few arguments!\n");
        return 1;
    }

    /* Processing arguments of command line */
    ROWS = atoi(argv[1]);
    COLUMNS = atoi(argv[2]);
    
    if (argc == 4) {
        fp = fopen(argv[3], "wb");
    } else {
        fp = fopen("data.txt", "wb");
    }

    /* Seting seed for a new sequence of pseudo-random integers */
    srand(time(NULL));
    fprintf(fp, "%d %d\n", ROWS, COLUMNS);
    /* Matrix generating */
    for (i = 0; i < ROWS; ++i) {
        for (j = 0; j < COLUMNS; ++j) {
            /* Sign of random number */
            if (rand() % 10 + 1 > 5) {
                fprintf(fp, "%f ", rand() / (double) RAND_MAX * 10);
            } else {
                fprintf(fp, "%f ", 0.0 - rand() / (double) RAND_MAX * 10);
            }
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
        
    return 0;
}
