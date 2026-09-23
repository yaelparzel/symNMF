#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

/*
 * Allocate a rows x cols matrix with every entry set to zero.
 * On failure releases whatever was already allocated and returns NULL.
 */
double **mat_alloc(int rows, int cols)
{
    double **m;
    int i;

    if (rows <= 0 || cols <= 0) {
        return NULL;
    }
    m = (double **)calloc((size_t)rows, sizeof(double *));
    if (m == NULL) {
        return NULL;
    }
    for (i = 0; i < rows; i++) {
        m[i] = (double *)calloc((size_t)cols, sizeof(double));
        if (m[i] == NULL) {
            mat_free(m, i);
            return NULL;
        }
    }
    return m;
}

/*
 * Return a deep copy of the (rows x cols) matrix m.
 * Returns NULL on allocation failure.
 */
double **mat_copy(double **m, int rows, int cols)
{
    double **c;
    int i, j;

    c = mat_alloc(rows, cols);
    if (c == NULL) {
        return NULL;
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            c[i][j] = m[i][j];
        }
    }
    return c;
}

/*
 * Compute C = A (ar x ac) * B (ac x bc) and return new ar x bc matrix.
 * Returns NULL on allocation failure; caller frees result with mat_free.
 * Assumes inner dimensions (ac) match; no runtime check performed.
 */
double **mat_mult(double **a, int ar, int ac, double **b, int bc)
{
    double **c;
    double sum;
    int i, j, t;

    c = mat_alloc(ar, bc);
    if (c == NULL) {
        return NULL;
    }
    for (i = 0; i < ar; i++) {
        for (j = 0; j < bc; j++) {
            sum = 0.0;
            for (t = 0; t < ac; t++) {
                sum += a[i][t] * b[t][j];
            }
            c[i][j] = sum;
        }
    }
    return c;
}

/*
 * Transpose of the (rows x cols) matrix m.
 * Returns a newly allocated (cols x rows) matrix, or NULL on failure.
 */
double **mat_transpose(double **m, int rows, int cols)
{
    double **t;
    int i, j;

    t = mat_alloc(cols, rows);
    if (t == NULL) {
        return NULL;
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            t[j][i] = m[i][j];
        }
    }
    return t;
}

/*
 * Squared Frobenius norm of elementwise (a - b) for rows x cols.
 * Returns sum of squared elementwise differences: (a[i][j]-b[i][j])^2.
 */
double mat_sq_frob_diff(double **a, double **b, int rows, int cols)
{
    double sum = 0.0;
    double diff;
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            diff = a[i][j] - b[i][j];
            sum += diff * diff;
        }
    }
    return sum;
}

/*
 * Squared Euclidean distance between vectors p and q, each of length d.
 * Returns sum over i of (p[i] - q[i])^2; no square root is taken.
 */
double sq_euclidean(const double *p, const double *q, int d)
{
    double sum = 0.0;
    double diff;
    int i;

    for (i = 0; i < d; i++) {
        diff = p[i] - q[i];
        sum += diff * diff;
    }
    return sum;
}

/*
 * Print the (rows x cols) matrix m: every value to four decimal places, values within a row separated by commas, one row per line.
 */
void mat_print(double **m, int rows, int cols)
{
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (j > 0) {
                printf(",");
            }
            printf("%.4f", m[i][j]);
        }
        printf("\n");
    }
}

/*
 * Free a matrix previously returned by mat_alloc; safe to call with NULL.
 * Frees the first row buffers, then the row pointer array.
 */
void mat_free(double **m, int rows)
{
    int i;

    if (m == NULL) {
        return;
    }
    for (i = 0; i < rows; i++) {
        free(m[i]);
    }
    free(m);
}