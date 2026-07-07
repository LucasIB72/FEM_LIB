#ifndef SOLVER_H
#define SOLVER_H

/* ------------------------------------------------------------------
   Dense solver (LAPACK dsysv-like interface)
   A[0..n-1][0..n-1] in column-major — converte internamente se
   necess�rio. Resolve A * X = B. A � sobrescrita.
   ------------------------------------------------------------------ */
void dsolve_sym(int n, int nrhs, double* A, double* B);

/* ------------------------------------------------------------------
   Solver esparso CSR via Gradiente Conjugado (CG).
   Resolve K * u = rhs. Retorna n�mero de itera��es ou -1 se n�o
   convergir.
   ------------------------------------------------------------------ */
int solve_csr_cg(int n, const int* rowIndex, const int* columns,
    const double* values, const double* rhs, double* u,
    int max_iter, double tol);

/* ------------------------------------------------------------------
   PARDISO stub (quando MKL estiver disponivel, substituir
   a implementacao).
   pt, iparm, mtype etc seguem a conven��o PARDISO 5.0.
   ------------------------------------------------------------------ */
int solve_csr_pardiso(int n, const int* rowIndex, const int* columns,
    const double* values, const double* rhs, double* u,
    int nrhs);

#endif
