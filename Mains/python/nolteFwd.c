#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <arrayobject.h>

#include <unistd.h>
#include <float.h>
#include <strings.h>
#include <fcntl.h>
#include <model.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include "samutil.h"

#include <gsl/gsl_matrix.h>

static char *DSName;            /* there can be only one */
static HeaderInfo Header;
static ChannelInfo *Channel;
static EpochInfo *Epoch;
static HULL Hull;
static COEFFS *Coeffs;
static int M;
static gsl_matrix *L;

/* Python wrapper code. */

static char Doc_dsopen[] = "dsopen(filename)\n\
Open a CTF dataset.\n";

static PyObject *dsopen_wrap(PyObject *self, PyObject *args)
{
    char *dsname, *s;
    int i;

    if (!PyArg_ParseTuple(args, "s", &dsname)) {
        return NULL;
    }

    DSName = (char *)malloc(strlen(dsname) + 1);
    strcpy(DSName, dsname);

    fprintf(stderr, "opening %s\n", DSName);

    GetDsInfo(DSName, &Header, &Channel, &Epoch, NULL, TRUE);

    /* read the hull.shape file */

    i = strlen(Header.DsPath) + 2;
    i += strlen(Header.SetName) + 15;
    s = (char *)malloc(i);
    strcpy(s, Header.DsPath);
    strcat(s, "/");
    strcat(s, Header.SetName);
    strcat(s, ".ds/hull.shape");
    GetHull(s, &Hull);
    free(s);

    ORDER = 16;
    Coeffs = new(COEFFS);
    ECDIntPnt(&Header, Channel, &Hull, Coeffs);

    /* allocate the leadfield array */

    M = Header.NumPri;
    L = gsl_matrix_alloc(M, 3);

    Py_INCREF(Py_None);
    return Py_None;
}

static char Doc_doFwd[] = "doFwd(pos, ori)\n\
Return the lead field for a dipole at pos with orientation ori.\n\
pos and ori are three-tuples specifying position in cm and orientation,\n\
respectively.\n";

static PyObject *doFwd_wrap(PyObject *self, PyObject *args)
{
    int i, j;
    double pos[3], ori[3], d, *b;
    npy_intp dim[1];
    PyArrayObject *r;

    if (!PyArg_ParseTuple(args, "(ddd)(ddd)", pos, pos + 1, pos + 2, ori, ori + 1, ori + 2)) {
        return NULL;
    }

    /* Convert to meters. */

    for (j = 0; j < 3; j++) {
        pos[j] *= .01;
    }

    /* Normalize the passed orientation vector. */

    d = 0.;
    for (j = 0; j < 3; j++) {
        d += ori[j] * ori[j];
    }
    d = sqrt(d);
    for (j = 0; j < 3; j++) {
        ori[j] /= d;
    }

    /* Calculate the lead field, L. */

    ECDLeadField0(pos, L, &Header, Channel, &Hull, Coeffs);

    /* Dot the lead field with the orientation to get the result. */

    dim[0] = M;
    r = (PyArrayObject *)PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    if (r == NULL) {
        return NULL;
    }
    b = (double *)PyArray_DATA(r);

    for (i = 0; i < M; i++) {
        d = 0.;
        for (j = 0; j < 3; j++) {
            d += ori[j] * gsl_matrix_get(L, i, j);
        }
        b[i] = d;
    }

    return PyArray_Return(r);
}

static char Doc_doFwdMoment[] = "doFwdMoment(pos, Cinv)\n\
Calculate the orientation using the data's inverse covariance matrix.\n\
Return the lead field for a dipole at pos with the calculated moment.\n\
pos is a three-tuple specifying position in cm, and Cinv is an MxM matrix.\n";

static PyObject *doFwdMoment_wrap(PyObject *self, PyObject *args)
{
    int i, j;
    double pos[3], ori[3], d, *b;
    npy_intp dim[1];
    PyObject *o;
    PyArrayObject *a, *r;
    gsl_matrix *Cinv;

    if (!PyArg_ParseTuple(args, "(ddd)O", pos, pos + 1, pos + 2, &o)) {
        return NULL;
    }

    /* Convert to meters. */

    for (j = 0; j < 3; j++) {
        pos[j] *= .01;
    }

    /* Get the inverse covariance matrix. */

    a = (PyArrayObject *)PyArray_ContiguousFromAny(o, NPY_DOUBLE, 2, 2);
    if (a == NULL) {
        return NULL;
    }
    i = PyArray_DIM(a, 0);
    j = PyArray_DIM(a, 1);
    if (i != j || i != M) {
        PyErr_SetString(PyExc_ValueError, "Cinv must be MxM!");
        Py_DECREF(a);
        return NULL;
    }
    Cinv = gsl_matrix_alloc(M, M);
    b = (double *)PyArray_DATA(a);
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            gsl_matrix_set(Cinv, i, j, *b++);
        }
    }
    Py_DECREF(a);

    /* Calculate the lead field, L. */

    ECDLeadField0(pos, L, &Header, Channel, &Hull, Coeffs);

    /* Solve generalized eigensystem for moment vector. */

    SolveMoment(Cinv, L, ori, &d);

    /* Dot the lead field with the orientation to get the result. */

    dim[0] = M;
    r = (PyArrayObject *)PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    if (r == NULL) {
        gsl_matrix_free(Cinv);
        return NULL;
    }
    b = (double *)PyArray_DATA(r);

    for (i = 0; i < M; i++) {
        d = 0.;
        for (j = 0; j < 3; j++) {
            d += ori[j] * gsl_matrix_get(L, i, j);
        }
        b[i] = d;
    }

    gsl_matrix_free(Cinv);
    return PyArray_Return(r);
}

static char Doc_nolteFwd[] = "Compute MEG forward solutions using the Nolte model.";

static PyMethodDef Methods[] = {
    { "dsopen", dsopen_wrap, METH_VARARGS, Doc_dsopen },
    { "doFwd", doFwd_wrap, METH_VARARGS, Doc_doFwd },
    { "doFwdMoment", doFwdMoment_wrap, METH_VARARGS, Doc_doFwdMoment },
    { NULL, NULL, 0, NULL }
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "nolteFwd",
    Doc_nolteFwd,
    -1,
    Methods,
    NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit_nolteFwd(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }

    import_array();

    return m;
}

#else

PyMODINIT_FUNC initnolteFwd()
{
    Py_InitModule3("nolteFwd", Methods, Doc_nolteFwd);
    import_array();
}

#endif
