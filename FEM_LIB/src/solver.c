#include "solver.h"
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

/* ==================================================================
   Solver denso: LAPACK (Intel MKL) via LAPACKE_dsysv
   ================================================================== */

void dsolve_sym(int n, int nrhs, double* A, double* B)
{
    lapack_int* ipiv = malloc(n * sizeof(lapack_int));
    if (!ipiv) return;

    /* A e B estao armazenadas em row-major (mesma convencao de
       assembly.c: K_global[i * n + j]); LAPACKE cuida da conversao
       interna para a ordem que a rotina Fortran do LAPACK espera. */
    lapack_int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', n, nrhs, A, n, ipiv, B, nrhs);

    if (info != 0)
        printf("LAPACK dsysv: falha na solucao (info = %d)\n", (int)info);

    free(ipiv);
}

/* ==================================================================
   Solver esparso: Intel MKL PARDISO (direto, paralelo)
   ================================================================== */

int solve_csr_pardiso(int n, const int* rowIndex, const int* columns,
    const double* values, const double* rhs, double* u, int nrhs)
{
    /* PARDISO exige indices em MKL_INT; converte os vetores CSR
       (rowIndex/columns), que na montagem sao "int" simples. */
    MKL_INT* ia = malloc((n + 1) * sizeof(MKL_INT));
    MKL_INT* ja = malloc(rowIndex[n] * sizeof(MKL_INT));
    if (!ia || !ja) { free(ia); free(ja); return -1; }

    for (int i = 0; i <= n; i++) ia[i] = rowIndex[i];
    for (int i = 0; i < rowIndex[n]; i++) ja[i] = columns[i];

    void* pt[64] = { 0 };
    MKL_INT iparm[64] = { 0 };
    MKL_INT mtype = 11;   /* matriz real, nao-simetrica (armazenamento CSR completo, ambos os triangulos) */
    MKL_INT maxfct = 1, mnum = 1, msglvl = 0, error = 0;
    MKL_INT n_mkl = n, nrhs_mkl = nrhs, idum = 0;

    pardisoinit(pt, &mtype, iparm);
    iparm[34] = 1;   /* indexacao 0-based (estilo C) em ia/ja */
    iparm[26] = 1;   /* verifica a estrutura da matriz CSR (ordem das colunas, duplicatas, etc.) antes de resolver */

    MKL_INT phase = 13; /* analise + fatoracao numerica + solucao + refinamento iterativo */
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n_mkl, (void*)values, ia, ja,
        &idum, &nrhs_mkl, iparm, &msglvl, (void*)rhs, u, &error);

    /* Libera as estruturas internas do PARDISO (usa uma variavel separada
       para nao sobrescrever o resultado da fase de solucao acima). */
    MKL_INT release_error = 0;
    phase = -1;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n_mkl, (void*)values, ia, ja,
        &idum, &nrhs_mkl, iparm, &msglvl, (void*)rhs, u, &release_error);

    free(ia);
    free(ja);

    return (int)error;
}
