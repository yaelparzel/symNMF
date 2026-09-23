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
 * Measure the layout of the input stream: *n receives the number of lines
 * holding data and *d the number of comma separated values on the first
 * such line. Blank lines are ignored and both LF and CRLF endings are
 * accepted, with or without a final newline.
 * The stream is rewound before returning. Returns 0 on success, 1 if the
 * file holds no data at all.
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
 * Consume the separator that follows a value: any run of whitespace, plus
 * one comma if one is present. Anything else is pushed back for the next
 * read, so the stream is left on the first character of the next value.
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
 * Read n points of d coordinates each from fp into a new matrix. Values are
 * read with fscanf, which skips leading whitespace, and the separator that
 * follows each value is consumed explicitly.
 * Returns NULL on allocation failure or on a value that cannot be parsed.
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
 * Read the data points of the named file, setting *n to the number of
 * points and *d to their dimension.
 * Returns NULL if the file cannot be opened, contains no data, or contains
 * something that is not a number.
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
 * Similarity matrix A (n x n) of the n data points of x, each of length d:
 * a_ij = exp(-||x_i - x_j||^2 / 2) for i != j, and a_ii = 0.
 * Returns a newly allocated matrix, or NULL on allocation failure.
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
 * Allocate and return the vector of row sums of the n x n matrix m, that
 * is the degree of every vertex of the similarity graph.
 * Returns NULL on allocation failure.
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
 * Diagonal degree matrix D (n x n) of the n data points of x, each of
 * length d. Entry (i,i) is the degree of point i and every off-diagonal
 * entry is zero.
 * Returns a newly allocated matrix, or NULL on failure.
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
 * Normalised similarity matrix W = D^(-1/2) A D^(-1/2) (n x n) for the n
 * data points of x, each of length d. Entry by entry this is
 * w_ij = a_ij / sqrt(deg_i * deg_j).
 * Returns a newly allocated matrix, or NULL on failure.
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
 * Build H(t+1) entrywise from H, the product W*H and the product H*(H^T H):
 * next_ij = h_ij * (1 - beta + beta * wh_ij / hhth_ij).
 * A zero denominator can only occur where h_ij is itself zero, so the entry
 * is set to zero rather than dividing.
 * Returns a newly allocated n x k matrix, or NULL on allocation failure.
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
 * Perform one multiplicative update of h against w.
 * Returns a newly allocated n x k matrix holding H(t+1), or NULL if any
 * intermediate allocation fails. All intermediates are released here.
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
 * Optimise H (n x k) against W (n x n), stopping once SYMNMF_MAX_ITER
 * updates have run or the squared Frobenius norm of the change between two
 * consecutive iterates drops below SYMNMF_EPS.
 * h is not modified; the final H is returned as a new matrix, NULL on failure.
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
 * Dispatch to the routine named by goal, which must be one of sym, ddg or
 * norm. The result is the corresponding n x n matrix.
 * Returns NULL on an unrecognised goal or on failure.
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
 * Entry point of the standalone program, invoked as
 * ./symnmf <goal> <file_name> with goal one of sym, ddg or norm.
 * Prints the requested matrix, or the error message on any failure.
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