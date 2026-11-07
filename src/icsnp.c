/*
 * Compiled routines of ICSNP, called through .Call.
 *
 * All matrices are passed and returned in R's column-major order, so the
 * element (i, e) of an n x k matrix is x[i + e * n].
 */

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <limits.h>

/* how often (in outer-loop iterations) to check for a user interrupt */
#define ICSNP_INTERRUPT_EVERY 64

static void check_real_matrix(SEXP X)
{
    if (!isReal(X) || !isMatrix(X))
        error("'X' must be a numeric matrix");
}

/*
 * Pairwise operations on the rows of X: for all i < j the row
 * x_i - x_j, x_i + x_j or x_i * x_j (elementwise). The pairs are ordered
 * (1,2), (1,3), ..., (1,n), (2,3), ..., (n-1,n).
 */

enum pair_op { PAIR_DIFF, PAIR_SUM, PAIR_PROD };

static SEXP pair_op(SEXP X, enum pair_op op)
{
    check_real_matrix(X);
    R_xlen_t n = nrows(X), k = ncols(X);
    R_xlen_t np = n > 1 ? n * (n - 1) / 2 : 0;
    if (np > INT_MAX)
        error("too many pairs: n * (n - 1) / 2 must not exceed %d", INT_MAX);
    const double *x = REAL(X);

    SEXP res = PROTECT(allocMatrix(REALSXP, (int) np, (int) k));
    double *r = REAL(res);

    for (R_xlen_t e = 0; e < k; e++) {
        const double *xe = x + e * n;
        double *re = r + e * np;
        R_xlen_t row = 0;
        for (R_xlen_t i = 0; i < n - 1; i++) {
            if (i % ICSNP_INTERRUPT_EVERY == 0) R_CheckUserInterrupt();
            const double xi = xe[i];
            switch (op) {
            case PAIR_DIFF:
                for (R_xlen_t j = i + 1; j < n; j++) re[row++] = xi - xe[j];
                break;
            case PAIR_SUM:
                for (R_xlen_t j = i + 1; j < n; j++) re[row++] = xi + xe[j];
                break;
            case PAIR_PROD:
                for (R_xlen_t j = i + 1; j < n; j++) re[row++] = xi * xe[j];
                break;
            }
        }
    }
    UNPROTECT(1);
    return res;
}

SEXP C_pairdiff(SEXP X) { return pair_op(X, PAIR_DIFF); }
SEXP C_pairsum(SEXP X)  { return pair_op(X, PAIR_SUM); }
SEXP C_pairprod(SEXP X) { return pair_op(X, PAIR_PROD); }

/*
 * Symmetrized shape and scatter matrices.
 *
 * Both Duembgen's shape matrix and the symmetrized Huber scatter matrix
 * are fixed points of a map of the form
 *
 *     V -> average over all pairs i < j of  w(r_ij) d_ij d_ij',
 *
 * with d_ij = x_i - x_j and r_ij = d_ij' V^{-1} d_ij, where w(r) = 1 / r
 * for Duembgen's shape matrix (followed by a standardization to
 * determinant one) and w(r) = 1 / sig if r <= c and (c / r) / sig otherwise
 * for the symmetrized Huber scatter matrix.
 *
 * One step is computed with the Cholesky factor V = L L': the data are
 * transformed once to z_i = L^{-1} x_i, so that with e_ij = z_i - z_j
 * r_ij = |e_ij|^2 and d_ij d_ij' = L e_ij e_ij' L', and the new matrix is
 * L S L' with S the average of w(r_ij) e_ij e_ij'.
 */

enum weight_type { WEIGHT_SIGN, WEIGHT_HUBER };

/* workspace, allocated once per call with R_alloc and freed by R */
typedef struct {
    int n, k;
    double *zt;  /* transformed data, k x n (one observation per column) */
    double *L;   /* Cholesky factor, k x k */
    double *S;   /* sum of weighted outer products, k x k */
    double *T;   /* k x k */
    double *e;   /* one pairwise difference, length k */
} workspace;

static workspace make_workspace(int n, int k)
{
    workspace w;
    w.n = n;
    w.k = k;
    w.zt = (double *) R_alloc((size_t) n * k, sizeof(double));
    w.L = (double *) R_alloc((size_t) k * k, sizeof(double));
    w.S = (double *) R_alloc((size_t) k * k, sizeof(double));
    w.T = (double *) R_alloc((size_t) k * k, sizeof(double));
    w.e = (double *) R_alloc(k, sizeof(double));
    return w;
}

/* copy the lower triangle of a k x k matrix to the upper one */
static void symmetrize(double *y, int k)
{
    for (int m = 0; m < k; m++)
        for (int q = m + 1; q < k; q++)
            y[m + q * k] = y[q + m * k];
}

/* lower Cholesky factor L of the symmetric k x k matrix a (the upper part
   of L is set to zero); returns 0 if a is not positive definite */
static int chol_lower(const double *a, double *L, int k)
{
    for (int m = 0; m < k * k; m++) L[m] = 0.0;
    for (int j = 0; j < k; j++) {
        double s = a[j + j * k];
        for (int q = 0; q < j; q++) s -= L[j + q * k] * L[j + q * k];
        if (!(s > 0.0)) return 0;
        const double ljj = sqrt(s);
        L[j + j * k] = ljj;
        for (int i = j + 1; i < k; i++) {
            double t = a[i + j * k];
            for (int q = 0; q < j; q++) t -= L[i + q * k] * L[j + q * k];
            L[i + j * k] = t / ljj;
        }
    }
    return 1;
}

static double frobenius_dist(const double *a, const double *b, int len)
{
    double s = 0.0;
    for (int m = 0; m < len; m++) s += (a[m] - b[m]) * (a[m] - b[m]);
    return sqrt(s);
}

/*
 * One step of the map above for the n x k data x (column-major). V is the
 * current matrix, or NULL for the identity. The result is written to out,
 * and the number of zero pairwise differences, which are skipped for
 * WEIGHT_SIGN, is returned.
 */
static double one_step(const double *x, const double *V, enum weight_type type,
                       double c, double sig, workspace *w, double *out)
{
    const int n = w->n, k = w->k;
    double *zt = w->zt, *L = w->L, *S = w->S, *T = w->T, *e = w->e;

    /* z_i = L^{-1} x_i by forward substitution */
    if (V != NULL) {
        if (!chol_lower(V, L, k))
            error("the current scatter matrix is not positive definite");
        for (int i = 0; i < n; i++) {
            double *zi = zt + (R_xlen_t) i * k;
            for (int m = 0; m < k; m++) {
                double s = x[i + (R_xlen_t) m * n];
                for (int q = 0; q < m; q++) s -= L[m + q * k] * zi[q];
                zi[m] = s / L[m + m * k];
            }
        }
    } else {
        for (int i = 0; i < n; i++)
            for (int m = 0; m < k; m++)
                zt[m + (R_xlen_t) i * k] = x[i + (R_xlen_t) m * n];
    }

    /* S = sum of w(r) e e' over the pairs, accumulated in the lower triangle */
    for (int m = 0; m < k * k; m++) S[m] = 0.0;
    double n_zero = 0.0;
    for (int i = 0; i < n - 1; i++) {
        if (i % ICSNP_INTERRUPT_EVERY == 0) R_CheckUserInterrupt();
        const double *zi = zt + (R_xlen_t) i * k;
        for (int j = i + 1; j < n; j++) {
            const double *zj = zt + (R_xlen_t) j * k;
            double r = 0.0;
            for (int m = 0; m < k; m++) {
                e[m] = zi[m] - zj[m];
                r += e[m] * e[m];
            }
            double wt;
            if (type == WEIGHT_SIGN) {
                /* a zero difference (tied observations) has no sign */
                if (r == 0.0) {
                    n_zero += 1.0;
                    continue;
                }
                wt = 1.0 / r;
            } else {
                wt = (r <= c) ? 1.0 / sig : (c / r) / sig;
            }
            for (int q = 0; q < k; q++) {
                const double weq = wt * e[q];
                for (int m = q; m < k; m++) S[m + q * k] += e[m] * weq;
            }
        }
    }
    symmetrize(S, k);
    const double n_pairs = (double) n * (n - 1) / 2.0;
    for (int m = 0; m < k * k; m++) S[m] /= n_pairs;

    if (V == NULL) {
        for (int m = 0; m < k * k; m++) out[m] = S[m];
        return n_zero;
    }
    /* out = L S L', using that L is lower triangular */
    for (int b = 0; b < k; b++)
        for (int a = 0; a < k; a++) {
            double s = 0.0;
            for (int q = 0; q <= a; q++) s += L[a + q * k] * S[q + b * k];
            T[a + b * k] = s;
        }
    for (int b = 0; b < k; b++)
        for (int a = b; a < k; a++) {
            double s = 0.0;
            for (int q = 0; q <= b; q++) s += T[a + q * k] * L[b + q * k];
            out[a + b * k] = s;
        }
    symmetrize(out, k);
    return n_zero;
}

/* scale the k x k matrix V to determinant one */
static void standardize_det(double *V, workspace *w)
{
    const int k = w->k;
    if (!chol_lower(V, w->L, k))
        error("the shape matrix is not positive definite");
    double logdet = 0.0;
    for (int m = 0; m < k; m++) logdet += 2.0 * log(w->L[m + m * k]);
    const double scale = exp(-logdet / k);
    for (int m = 0; m < k * k; m++) V[m] *= scale;
}

static void check_square(SEXP V, int k, const char *name)
{
    if (!isReal(V) || !isMatrix(V) || nrows(V) != k || ncols(V) != k)
        error("'%s' must be a numeric %d x %d matrix", name, k, k);
}

/*
 * Average over all pairs i < j of the outer products of the spatial signs
 * of x_i - x_j. Pairs with a zero difference (tied observations) are
 * skipped, but counted in the number of pairs; their number is returned as
 * the attribute "n.zero".
 */
SEXP C_sum_diff_sign_outers(SEXP X)
{
    check_real_matrix(X);
    const int n = nrows(X), k = ncols(X);
    workspace w = make_workspace(n, k);
    SEXP res = PROTECT(allocMatrix(REALSXP, k, k));
    double n_zero = one_step(REAL(X), NULL, WEIGHT_SIGN, 0.0, 0.0, &w, REAL(res));
    SEXP nz = PROTECT(ScalarReal(n_zero));
    setAttrib(res, install("n.zero"), nz);
    UNPROTECT(2);
    return res;
}

/* one step of the symmetrized Huber iteration from the scatter matrix V */
SEXP C_symm_huber_step(SEXP X, SEXP V, SEXP cs, SEXP sigs)
{
    check_real_matrix(X);
    const int n = nrows(X), k = ncols(X);
    check_square(V, k, "V");
    workspace w = make_workspace(n, k);
    SEXP res = PROTECT(allocMatrix(REALSXP, k, k));
    one_step(REAL(X), REAL(V), WEIGHT_HUBER, asReal(cs), asReal(sigs), &w, REAL(res));
    UNPROTECT(1);
    return res;
}

/*
 * Duembgen's shape matrix by the fixed-point iteration, started from V0
 * (which should have determinant one). If steps is finite, exactly that
 * many steps are taken; otherwise the iteration stops when the Frobenius
 * distance between two iterates is below eps. The attribute "iter" gives
 * the number of steps.
 */
SEXP C_duembgen_shape(SEXP X, SEXP V0, SEXP Reps, SEXP Rmaxiter, SEXP Rsteps)
{
    check_real_matrix(X);
    const int n = nrows(X), k = ncols(X);
    check_square(V0, k, "init");
    const double eps = asReal(Reps), maxiter = asReal(Rmaxiter), steps = asReal(Rsteps);
    workspace w = make_workspace(n, k);
    double *Vcur = (double *) R_alloc((size_t) k * k, sizeof(double));
    double *Vnew = (double *) R_alloc((size_t) k * k, sizeof(double));
    for (int m = 0; m < k * k; m++) Vcur[m] = REAL(V0)[m];

    int iter = 0;
    while (TRUE) {
        if (iter >= steps) break;
        if (iter >= maxiter) error("maxiter reached");
        iter++;
        one_step(REAL(X), Vcur, WEIGHT_SIGN, 0.0, 0.0, &w, Vnew);
        standardize_det(Vnew, &w);
        const int converged = !R_FINITE(steps) && frobenius_dist(Vnew, Vcur, k * k) < eps;
        for (int m = 0; m < k * k; m++) Vcur[m] = Vnew[m];
        if (converged) break;
    }

    SEXP res = PROTECT(allocMatrix(REALSXP, k, k));
    for (int m = 0; m < k * k; m++) REAL(res)[m] = Vcur[m];
    SEXP it = PROTECT(ScalarInteger(iter));
    setAttrib(res, install("iter"), it);
    UNPROTECT(2);
    return res;
}

/*
 * The symmetrized Huber scatter matrix by the fixed-point iteration,
 * started from V0, until the Frobenius distance between two iterates is
 * below eps. The attribute "iter" gives the number of steps.
 */
SEXP C_symm_huber(SEXP X, SEXP V0, SEXP cs, SEXP sigs, SEXP Reps, SEXP Rmaxiter)
{
    check_real_matrix(X);
    const int n = nrows(X), k = ncols(X);
    check_square(V0, k, "init");
    const double csq = asReal(cs), sigsq = asReal(sigs);
    const double eps = asReal(Reps), maxiter = asReal(Rmaxiter);
    workspace w = make_workspace(n, k);
    double *Vcur = (double *) R_alloc((size_t) k * k, sizeof(double));
    double *Vnew = (double *) R_alloc((size_t) k * k, sizeof(double));
    for (int m = 0; m < k * k; m++) Vcur[m] = REAL(V0)[m];

    int iter = 0;
    double differ = R_PosInf;
    while (differ > eps) {
        if (iter >= maxiter) error("maxiter reached without convergence");
        iter++;
        one_step(REAL(X), Vcur, WEIGHT_HUBER, csq, sigsq, &w, Vnew);
        differ = frobenius_dist(Vnew, Vcur, k * k);
        for (int m = 0; m < k * k; m++) Vcur[m] = Vnew[m];
    }

    SEXP res = PROTECT(allocMatrix(REALSXP, k, k));
    for (int m = 0; m < k * k; m++) REAL(res)[m] = Vcur[m];
    SEXP it = PROTECT(ScalarInteger(iter));
    setAttrib(res, install("iter"), it);
    UNPROTECT(2);
    return res;
}
