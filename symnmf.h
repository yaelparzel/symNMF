#ifndef SYMNMF_H_
#define SYMNMF_H_

/*
 * X is the (n x d) matrix of data points.
 * all returned matrices are allocated and must be released with mat_free. NULL is returned on failure.
 */

/* Similarity matrix A (n x n): a_ij = exp(-||x_i - x_j||^2 / 2), a_ii = 0. */
double **sym(double **x, int n, int d);

/* Diagonal degree matrix D (n x n), row sums of A on the diagonal. */
double **ddg(double **x, int n, int d);

/* Normalised similarity W = D^(-1/2) A D^(-1/2) (n x n). */
double **norm(double **x, int n, int d);

/*
 * Optimise H (n x k) against W (n x n) with the multiplicative update rule
 * until SYMNMF_MAX_ITER is reached or ||H(t+1) - H(t)||_F^2 < SYMNMF_EPS.
 * h is not modified, the final H is returned as a new matrix.
 */
double **symnmf(double **h, double **w, int n, int k);

#endif