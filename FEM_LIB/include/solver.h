#ifndef SOLVER_H
#define SOLVER_H

/* ------------------------------------------------------------------
   Solver denso via LAPACK (Intel MKL, LAPACKE_dsysv).
   Resolve A * X = B para uma matriz simetrica geral (nao
   necessariamente definida positiva). A e B (armazenadas em
   row-major, mesma convencao usada em assembly.c) sao sobrescritas
   com a fatoracao e a solucao, respectivamente.
   ------------------------------------------------------------------ */
void dsolve_sym(int n, int nrhs, double* A, double* B);

/* ------------------------------------------------------------------
   Solver esparso via Intel MKL PARDISO (direto, paralelo).
   Resolve K * u = rhs para uma matriz CSR (rowIndex/columns/values,
   indexacao 0-based). Retorna 0 em sucesso, ou o codigo de erro do
   PARDISO (negativo) em caso de falha.
   ------------------------------------------------------------------ */
int solve_csr_pardiso(int n, const int* rowIndex, const int* columns,
    const double* values, const double* rhs, double* u, int nrhs);

#endif
