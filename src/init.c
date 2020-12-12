#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following name(s) appear with different usages
  e.g., with different numbers of arguments:

    find_best_candidate

  This needs to be resolved in the tables and any declarations.
*/

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void find_best_candidate(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_I_vec(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"find_best_candidate", (DL_FUNC) &find_best_candidate, 26},
    {"get_I_vec",           (DL_FUNC) &get_I_vec,            7},
    {NULL, NULL, 0}
};

void R_init_CoxBoost(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
