#include <stdio.h>
#include <stdlib.h>
#include <Python.h>

typedef struct {
	PyObject_HEAD
	FILE* f;
	int res;
	char name1[100];
	char seq[400];
	char name2[100];
	char qual[400];
} cFileGenState;

PyObject* cFileGen_next(PyObject *self) {
	cFileGenState *p = (cFileGenState *)self;
	p->res = fscanf(p->f, "%s\n%s\n%s\n%s\n", p->name1, p->seq, p->name2, p->qual);
	if (p->res != EOF) {
		PyObject *tmp = Py_BuildValue("(ss)", p->seq, p->qual);
		return tmp;
	} else {
		PyErr_SetNone(PyExc_StopIteration);
		return NULL;
	}
}

static PyTypeObject cFileGenStateType = {
	PyObject_HEAD_INIT(NULL)
	0,
	"cIO._cFileGen",                /*tp_name*/
	sizeof(cFileGenState),     /*tp_basicsize*/
	0,0,0,0,0,0,0,0,                         /*tp_itemsize*/
	0,0,0,0,0,0,0,0,                         /*tp_as_number*/
	Py_TPFLAGS_DEFAULT,
	0,0,0,0,0,           /* tp_doc */
	PyObject_SelfIter,  /* tp_iter: __iter__() method */
	//cFileGen_iter,  /* tp_iter: __iter__() method */
	cFileGen_next  /* tp_iternext: next() method */
};

static PyObject* cFileGen(PyObject *self, PyObject *args) {
	char *fName=NULL;
	if (!PyArg_ParseTuple(args, "s", &fName)) return NULL;
	cFileGenState *p = PyObject_New(cFileGenState, &cFileGenStateType);
	if(!p) return NULL;
	p->f = fopen(fName,"r");
	return (PyObject *)p;
}

static PyMethodDef cFileGenMethods[] = {
    {"cFileGen", cFileGen, METH_VARARGS, "Fastq gen in C."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcIO(void) {
	PyObject *m;
	cFileGenStateType.tp_new = PyType_GenericNew;
	if(PyType_Ready(&cFileGenStateType) < 0) return;
	m = Py_InitModule("cIO", cFileGenMethods);
	if(!m) return;
	Py_INCREF(&cFileGenStateType);
	PyModule_AddObject(m,"_cFileGen", (PyObject*)&cFileGenStateType);
}
