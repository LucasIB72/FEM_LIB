#include "solver.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ==================================================================
   Dense symmetric solver  (dsysv-like)
   ================================================================== */

static void lu_decompose(int n, double* A, int* pivot)
{
    for (int k = 0; k < n; k++)
    {
        int p = k;
        double max_val = fabs(A[k * n + k]);
        for (int i = k + 1; i < n; i++)
        {
            if (fabs(A[i * n + k]) > max_val)
            {
                max_val = fabs(A[i * n + k]);
                p = i;
            }
        }
        pivot[k] = p;

        if (p != k)
        {
            double* row_k = A + k * n;
            double* row_p = A + p * n;
            for (int j = 0; j < n; j++)
            {
                double tmp = row_k[j];
                row_k[j] = row_p[j];
                row_p[j] = tmp;
            }
        }

        double inv_pivot = 1.0 / A[k * n + k];
        for (int i = k + 1; i < n; i++)
        {
            A[i * n + k] *= inv_pivot;
            for (int j = k + 1; j < n; j++)
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
        }
    }
}

static void lu_solve(int n, const double* A, const int* pivot, double* b)
{
    for (int i = 0; i < n; i++)
    {
        int p = pivot[i];
        if (p != i) { double tmp = b[i]; b[i] = b[p]; b[p] = tmp; }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            b[i] -= A[i * n + j] * b[j];

    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < n; j++)
            b[i] -= A[i * n + j] * b[j];
        b[i] /= A[i * n + i];
    }
}

void dsolve_sym(int n, int nrhs, double* A, double* B)
{
    int* pivot = malloc(n * sizeof(int));
    if (!pivot) return;

    lu_decompose(n, A, pivot);

    for (int r = 0; r < nrhs; r++)
        lu_solve(n, A, pivot, B + r * n);

    free(pivot);
}

/* ==================================================================
   Sparse CG solver
   ================================================================== */

static void csr_matvec(int n, const int* rowIndex, const int* columns,
    const double* values, const double* x, double* y)
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int p = rowIndex[i]; p < rowIndex[i + 1]; p++)
            sum += values[p] * x[columns[p]];
        y[i] = sum;
    }
}

static double dot(int n, const double* a, const double* b)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

static void axpy(int n, double alpha, const double* x, double* y)
{
    for (int i = 0; i < n; i++) y[i] += alpha * x[i];
}

static void xpay(int n, double alpha, const double* x, double* y)
{
    for (int i = 0; i < n; i++) y[i] = x[i] + alpha * y[i];
}

int solve_csr_cg(int n, const int* rowIndex, const int* columns,
    const double* values, const double* rhs, double* u,
    int max_iter, double tol)
{
    double* r = calloc(n, sizeof(double));
    double* p = calloc(n, sizeof(double));
    double* Ap = calloc(n, sizeof(double));

    if (!r || !p || !Ap) { free(r); free(p); free(Ap); return -1; }

    csr_matvec(n, rowIndex, columns, values, u, r);
    for (int i = 0; i < n; i++) r[i] = rhs[i] - r[i];

    memcpy(p, r, n * sizeof(double));

    double rr = dot(n, r, r);
    double tol_sq = tol * tol;

    int iter;
    for (iter = 0; iter < max_iter; iter++)
    {
        if (rr < tol_sq) break;

        csr_matvec(n, rowIndex, columns, values, p, Ap);

        double pAp = dot(n, p, Ap);
        if (fabs(pAp) < 1e-30) break;

        double alpha = rr / pAp;
        axpy(n, alpha, p, u);
        axpy(n, -alpha, Ap, r);

        double rr_new = dot(n, r, r);
        double beta = rr_new / rr;
        xpay(n, beta, r, p);
        rr = rr_new;
    }

    free(r); free(p); free(Ap);
    return iter;
}

/* ==================================================================
   PARDISO stub — compila mas retorna erro.
   Quando MKL estiver dispon�vel, substituir este .c pelo c�digo
   real que chama pardiso() da Intel.
   ================================================================== */

int solve_csr_pardiso(int n, const int* rowIndex, const int* columns,
    const double* values, const double* rhs, double* u,
    int nrhs)
{
    (void)n; (void)rowIndex; (void)columns; (void)values;
    (void)rhs; (void)u; (void)nrhs;
    return -1; /* PARDISO n�o dispon�vel */
}
