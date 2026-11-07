#include <stdlib.h> /* for NULL */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_pairdiff(SEXP);
extern SEXP C_pairsum(SEXP);
extern SEXP C_pairprod(SEXP);
extern SEXP C_sum_diff_sign_outers(SEXP);
extern SEXP C_symm_huber_step(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_duembgen_shape(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_symm_huber(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_pairdiff",             (DL_FUNC) &C_pairdiff,             1},
    {"C_pairsum",              (DL_FUNC) &C_pairsum,              1},
    {"C_pairprod",             (DL_FUNC) &C_pairprod,             1},
    {"C_sum_diff_sign_outers", (DL_FUNC) &C_sum_diff_sign_outers, 1},
    {"C_symm_huber_step",      (DL_FUNC) &C_symm_huber_step,      4},
    {"C_duembgen_shape",       (DL_FUNC) &C_duembgen_shape,       5},
    {"C_symm_huber",           (DL_FUNC) &C_symm_huber,           6},
    {NULL, NULL, 0}
};

void R_init_ICSNP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
