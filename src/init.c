/*
 * tree/src/init.c Copyright (C) 2002-12 B. D. Ripley
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 or 3 of the License
 *  (at your option).
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 */

#include <stddef.h>
#include <R.h>
#include <Rinternals.h>

#include "tree.h"
#include <R_ext/Rdynload.h>
#include <Rversion.h>

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
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
