#include <R.h>
#include <Rinternals.h>

#include "tree.h"
#include <R_ext/Rdynload.h>

#define CDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CMethodDef CEntries[]  = {
    CDEF(BDRgrow1, 23),
    CDEF(VR_dev1, 12),
    CDEF(VR_dev2, 10),
    CDEF(VR_dev3, 10),
    CDEF(VR_prune2, 17),
    CDEF(VR_pred1, 11),
    CDEF(VR_pred2, 10),
    {NULL, NULL, 0}
};


void R_init_tree(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
