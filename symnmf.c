#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "symnmf.h"

#define ERROR_MESSAGE "An Error Has Occurred"
#define SYMNMF_BETA 0.5
#define SYMNMF_EPS 1e-4
#define SYMNMF_MAX_ITER 300

/*
 * Determines the matrix dimensions from the input file.
 * Input: file pointer, pointers to row and column counters.
 * Output: 0 on success, or 1 if the file contains no data.
 */
static int scan_dimensions(FILE *fp, int *n, int *d)
{
    int c;
    int on_line = 0;
    int first_line = 1;

    *n = 0;
    *d = 1;
    while ((c = fgetc(fp)) != EOF) {
        if (c == '\n') {
            if (on_line) {
                (*n)++;
                first_line = 0;
            }
            on_line = 0;
        } 
        else if (c != '\r') {
            on_line = 1;
            if (c == ',' && first_line) {
                /* Count commas only on the first row to determine column count. */
                (*d)++;
            }
        }
    }
    if (on_line) {
        (*n)++;
    }
    rewind(fp);
    return (*n == 0) ? 1 : 0;
}

/*
 * Advances past the separator after a value.
 * Input: file pointer positioned after a parsed number.
 * Output: Leaves the stream ready to read the next value.
 */
static void skip_separator(FILE *fp)
{
    int c;

    do {
        c = fgetc(fp);
    } while (c == ' ' || c == '\t' || c == '\r' || c == '\n');
    if (c != ',' && c != EOF) {
        ungetc(c, fp);
    }
}

/*
 * Reads the data points into a matrix.
 * Input: file pointer, number of points, and dimension.
 * Output: Newly allocated matrix of values, or NULL on failure.
 */
static double **read_values(FILE *fp, int n, int d)
{
    double **x;
    int i, j;

    x = mat_alloc(n, d);
    if (x == NULL) {
        return NULL;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            if (fscanf(fp, "%lf", &x[i][j]) != 1) {
                mat_free(x, n);
                return NULL;
            }
            skip_separator(fp);
        }
    }
    return x;
}

/*
 * Loads the input points from a file.
 * Input: file name and pointers to row and column counts.
 * Output: Matrix of points, or NULL if the file is invalid.
 */
static double **read_points(const char *filename, int *n, int *d)
{
    FILE *fp;
    double **x;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        return NULL;
    }
    if (scan_dimensions(fp, n, d) != 0) {
        fclose(fp);
        return NULL;
    }
    x = read_values(fp, *n, *d);
    fclose(fp);
    return x;
}

/*
 * Computes the similarity matrix for the given data points.
 * Input: point matrix, number of points, and dimension.
 * Output: New n x n similarity matrix, or NULL on failure.
 */
double **sym(double **x, int n, int d)
{
    double **a;
    double value;
    int i, j;

    a = mat_alloc(n, n);
    if (a == NULL) {
        return NULL;
    }
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            value = exp(-sq_euclidean(x[i], x[j], d) / 2.0);
            a[i][j] = value;
            a[j][i] = value;
        }
    }
    return a;
}

/*
 * Computes the degree of each node in a matrix.
 * Input: matrix and number of rows.
 * Output: Array of row sums, or NULL on failure.
 */
static double *row_sums(double **m, int n)
{
    double *deg;
    int i, j;

    deg = (double *)calloc((size_t)n, sizeof(double));
    if (deg == NULL) {
        return NULL;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            deg[i] += m[i][j];
        }
    }
    return deg;
}

/*
 * Builds the diagonal degree matrix for the data.
 * Input: point matrix, number of points, and dimension.
 * Output: New degree matrix, or NULL on failure.
 */
double **ddg(double **x, int n, int d)
{
    double **a;
    double **deg_mat;
    double *deg;
    int i;

    a = sym(x, n, d);
    if (a == NULL) {
        return NULL;
    }
    deg = row_sums(a, n);
    deg_mat = mat_alloc(n, n);
    mat_free(a, n);
    if (deg == NULL || deg_mat == NULL) {
        free(deg);
        mat_free(deg_mat, n);
        return NULL;
    }
    for (i = 0; i < n; i++) {
        deg_mat[i][i] = deg[i];
    }
    free(deg);
    return deg_mat;
}

/*
 * Normalizes the similarity matrix using the degree matrix.
 * Input: point matrix, number of points, and dimension.
 * Output: New normalized matrix, or NULL on failure.
 */
double **norm(double **x, int n, int d)
{
    double **a;
    double **w;
    double *deg;
    int i, j;

    a = sym(x, n, d);
    if (a == NULL) {
        return NULL;
    }
    deg = row_sums(a, n);
    w = mat_alloc(n, n);
    if (deg == NULL || w == NULL) {
        free(deg);
        mat_free(w, n);
        mat_free(a, n);
        return NULL;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (deg[i] > 0.0 && deg[j] > 0.0) {
                w[i][j] = a[i][j] / sqrt(deg[i] * deg[j]);
            } else {
                w[i][j] = 0.0;
            }
        }
    }
    free(deg);
    mat_free(a, n);
    return w;
}

/*
 * Updates the matrix H using the multiplicative formula.
 * Input: current H, W*H product, and H*(H^T H) product.
 * Output: New H matrix, or NULL on failure.
 */
static double **fill_update(double **h, double **wh, double **hhth,
                            int n, int k)
{
    double **next;
    int i, j;

    next = mat_alloc(n, k);
    if (next == NULL) {
        return NULL;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            if (hhth[i][j] > 0.0) {
                /* Avoid division by zero when the denominator is zero. */
                next[i][j] = h[i][j] * (1.0 - SYMNMF_BETA +
                    SYMNMF_BETA * wh[i][j] / hhth[i][j]);
            } else {
                next[i][j] = 0.0;
            }
        }
    }
    return next;
}

/*
 * Performs one optimization step for H.
 * Input: current H matrix and normalized similarity matrix W.
 * Output: Updated H matrix, or NULL on failure.
 */
static double **update_h(double **h, double **w, int n, int k)
{
    double **wh = NULL, **ht = NULL, **hth = NULL;
    double **hhth = NULL, **next = NULL;

    wh = mat_mult(w, n, n, h, k);
    ht = mat_transpose(h, n, k);
    if (wh != NULL && ht != NULL) {
        hth = mat_mult(ht, k, n, h, k);
    }
    if (hth != NULL) {
        hhth = mat_mult(h, n, k, hth, k);
    }
    if (hhth != NULL) {
        next = fill_update(h, wh, hhth, n, k);
    }
    mat_free(wh, n);
    mat_free(ht, k);
    mat_free(hth, k);
    mat_free(hhth, n);
    return next;
}

/*
 * Optimizes the factor matrix H until convergence.
 * Input: initial H matrix, similarity matrix W, and dimensions.
 * Output: Final H matrix, or NULL on failure.
 */
double **symnmf(double **h, double **w, int n, int k)
{
    double **cur;
    double **next;
    int iter;

    cur = mat_copy(h, n, k);
    if (cur == NULL) {
        return NULL;
    }
    for (iter = 0; iter < SYMNMF_MAX_ITER; iter++) {
        next = update_h(cur, w, n, k);
        if (next == NULL) {
            mat_free(cur, n);
            return NULL;
        }
        if (mat_sq_frob_diff(next, cur, n, k) < SYMNMF_EPS) {
            mat_free(cur, n);
            return next;
        }
        mat_free(cur, n);
        cur = next;
    }
    return cur;
}

/*
 * Dispatches to the requested matrix-building function.
 * Input: goal name, data matrix, and dimensions.
 * Output: Result matrix for the selected goal, or NULL on failure.
 */
static double **run_goal(const char *goal, double **x, int n, int d)
{
    if (strcmp(goal, "sym") == 0) {
        return sym(x, n, d);
    }
    if (strcmp(goal, "ddg") == 0) {
        return ddg(x, n, d);
    }
    if (strcmp(goal, "norm") == 0) {
        return norm(x, n, d);
    }
    return NULL;
}

/*
 * Runs the command-line program for the selected matrix goal.
 * Input: command-line arguments with goal and input file.
 * Output: Prints the requested matrix or the error message.
 */
int main(int argc, char **argv)
{
    double **x;
    double **result;
    int n, d;

    if (argc != 3) {
        printf("%s\n", ERROR_MESSAGE);
        return 1;
    }
    x = read_points(argv[2], &n, &d);
    if (x == NULL) {
        printf("%s\n", ERROR_MESSAGE);
        return 1;
    }
    result = run_goal(argv[1], x, n, d);
    mat_free(x, n);
    if (result == NULL) {
        printf("%s\n", ERROR_MESSAGE);
        return 1;
    }
    mat_print(result, n, n);
    mat_free(result, n);
    return 0;
}