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
 * Return an independent copy of the rows x cols matrix m.
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
 * Matrix product a x b, where a is ar x ac and b is ac x bc.
 * The caller is responsible for the inner dimensions agreeing.
 * Returns a newly allocated ar x bc matrix, or NULL on allocation failure.
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
 * Transpose of the rows x cols matrix m.
 * Returns a newly allocated cols x rows matrix, or NULL on failure.
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
 * Squared Frobenius norm of (a - b), both rows x cols:
 * the sum of the squares of all entrywise differences.
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
 * Print the rows x cols matrix m: every value to four decimal places,
 * values within a row separated by commas, one row per line.
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
 * Release a matrix allocated by mat_alloc, freeing its first `rows` rows.
 * Safe to call on NULL.
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