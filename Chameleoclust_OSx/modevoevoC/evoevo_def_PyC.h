#ifndef PYC
#define PYC
#include <Python.h>
/**
* \struct PyCapsule
* \brief Contains a void* pointer (c object pointed), a char* name and a void* context and the pointer detructor
*/
typedef struct {
	PyObject_HEAD
	void *pointer;
	const char *name;
	void *context;
	PyCapsule_Destructor destructor;
} PyCapsule;
#endif
