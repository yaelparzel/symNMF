#ifndef MATRIX_H_
#define MATRIX_H_

/*
 * Allocate a rows x cols matrix initialized by zeros.
 * Returns NULL on allocation failure or if one of the dimensions is not positive.
 */
double **mat_alloc(int rows, int cols);

/* Free a matrix previously returned by a mat_* function. */
void mat_free(double **m, int rows);

/* Return an independent copy of the (rows x cols) matrix m. NULL on failure. */
double **mat_copy(double **m, int rows, int cols);

/* Print m, %.4f, values comma separated, one row per line. */
void mat_print(double **m, int rows, int cols);

/* Return (a x b), where a is (ar x ac) and b is (ac x bc). NULL on failure. */
double **mat_mult(double **a, int ar, int ac, double **b, int bc);

/* Return the transpose of m (cols x rows). NULL on failure. */
double **mat_transpose(double **m, int rows, int cols);

/* Return the squared Frobenius norm of (a - b). */
double mat_sq_frob_diff(double **a, double **b, int rows, int cols);

/* Return the squared Euclidean distance between vectors p and q of length d. */
double sq_euclidean(const double *p, const double *q, int d);

#endif