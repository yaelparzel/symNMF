#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include "matrix.h"
#include "symnmf.h"

/* Helper to convert a Python list of lists to a C 2D double array */
static double **py_list_to_c_matrix(PyObject *py_list, int *rows, int *cols) {
    int i, j;
    PyObject *row, *item;
    double **mat;

    if (!PyList_Check(py_list)) return NULL;

    *rows = (int)PyList_Size(py_list);
    if (*rows == 0) return NULL;

    row = PyList_GetItem(py_list, 0);
    if (!PyList_Check(row)) return NULL;
    
    *cols = (int)PyList_Size(row);
    if (*cols == 0) return NULL;

    mat = (double **)malloc((*rows) * sizeof(double *));
    if (!mat) return NULL;

    for (i = 0; i < *rows; i++) {
        row = PyList_GetItem(py_list, i);
        if (!PyList_Check(row) || (int)PyList_Size(row) != *cols) {
            mat_free(mat, i);
            return NULL;
        }
        mat[i] = (double *)malloc((*cols) * sizeof(double));
        if (!mat[i]) {
            mat_free(mat, i);
            return NULL;
        }
        for (j = 0; j < *cols; j++) {
            item = PyList_GetItem(row, j);
            mat[i][j] = PyFloat_AsDouble(item);
        }
    }
    if (PyErr_Occurred()) {
        mat_free(mat, *rows);
        return NULL;
    }
    return mat;
}

/* Helper to convert a C 2D double array to a Python list of lists */
static PyObject *c_matrix_to_py_list(double **mat, int rows, int cols) {
    int i, j;
    PyObject *py_list = PyList_New(rows);
    PyObject *py_row;

    if (!py_list) return NULL;

    for (i = 0; i < rows; i++) {
        py_row = PyList_New(cols);
        if (!py_row) {
            Py_DECREF(py_list);
            return NULL;
        }
        for (j = 0; j < cols; j++) {
            PyList_SetItem(py_row, j, PyFloat_FromDouble(mat[i][j]));
        }
        PyList_SetItem(py_list, i, py_row);
    }
    return py_list;
}

/* Python wrapper for sym(x, n, d) */
static PyObject* py_sym(PyObject *self, PyObject *args) {
    PyObject *py_x;
    double **x, **result;
    int n, d;
    PyObject *py_result;

    if (!PyArg_ParseTuple(args, "O", &py_x)) {
        return NULL;
    }

    x = py_list_to_c_matrix(py_x, &n, &d);
    if (!x) return NULL;

    /* Call the core logic function */
    result = sym(x, n, d);
    if (!result) {
        mat_free(x, n);
        return NULL;
    }

    py_result = c_matrix_to_py_list(result, n, n);

    mat_free(x, n);
    mat_free(result, n);

    return py_result;
}

/* Python wrapper for ddg(x, n, d) */
static PyObject* py_ddg(PyObject *self, PyObject *args) {
    PyObject *py_x;
    double **x, **result;
    int n, d;
    PyObject *py_result;

    if (!PyArg_ParseTuple(args, "O", &py_x)) {
        return NULL;
    }

    x = py_list_to_c_matrix(py_x, &n, &d);
    if (!x) return NULL;

    /* Call the core logic function */
    result = ddg(x, n, d);
    if (!result) {
        mat_free(x, n);
        return NULL;
    }
    
    py_result = c_matrix_to_py_list(result, n, n);

    mat_free(x, n);
    mat_free(result, n);

    return py_result;
}

/* Python wrapper for norm(x, n, d) */
static PyObject* py_norm(PyObject *self, PyObject *args) {
    PyObject *py_x;
    double **x, **result;
    int n, d;
    PyObject *py_result;

    if (!PyArg_ParseTuple(args, "O", &py_x)) {
        return NULL;
    }

    x = py_list_to_c_matrix(py_x, &n, &d);
    if (!x) return NULL;

    /* Call the core logic function */
    result = norm(x, n, d);
    if (!result) {
        mat_free(x, n);
        return NULL;
    }
    
    py_result = c_matrix_to_py_list(result, n, n);

    mat_free(x, n);
    mat_free(result, n);

    return py_result;
}

/* Python wrapper for symnmf(h, w, n, k) */
static PyObject* py_symnmf(PyObject *self, PyObject *args) {
    PyObject *py_h, *py_w;
    double **h, **w, **result;
    int n, k, n_w, n_cols_w;
    PyObject *py_result;

    /* Expecting two arguments: initial H matrix and normalized W matrix */
    if (!PyArg_ParseTuple(args, "OO", &py_h, &py_w)) {
        return NULL;
    }

    h = py_list_to_c_matrix(py_h, &n, &k);
    if (!h) return NULL;

    w = py_list_to_c_matrix(py_w, &n_w, &n_cols_w);
    if (!w) {
        mat_free(h, n);
        return NULL;
    }

    if (n_w != n || n_cols_w != n) {
        mat_free(h, n);
        mat_free(w, n_w);
        return NULL;
    }

    /* Call the core logic function */
    result = symnmf(h, w, n, k);
    if (!result) {
        mat_free(h, n);
        mat_free(w, n_w);
        return NULL;
    }
    
    py_result = c_matrix_to_py_list(result, n, k);

    mat_free(h, n);
    mat_free(w, n);
    mat_free(result, n);

    return py_result;
}

/* Map Python function names to C wrapper functions */
static PyMethodDef symnmfMethods[] = {
    {"sym", py_sym, METH_VARARGS, "Calculate the similarity matrix A"},
    {"ddg", py_ddg, METH_VARARGS, "Calculate the diagonal degree matrix D"},
    {"norm", py_norm, METH_VARARGS, "Calculate the normalized similarity matrix W"},
    {"symnmf", py_symnmf, METH_VARARGS, "Optimize H using the symNMF algorithm"},
    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    "Python C API for SymNMF",
    -1,
    symnmfMethods
};

/* Module initialization */
PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    return PyModule_Create(&symnmfmodule);
}