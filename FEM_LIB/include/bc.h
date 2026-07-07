#ifndef BC_H
#define BC_H

typedef struct {
    int n_dofs;
    int* dof;
    double* value;
} DisplacementBC;

typedef struct {
    int n_dofs;
    int* dof;
    double* value;
} ForceBC;

void apply_displacement_bc_dense(double* K, double* F, int n, const DisplacementBC* bc);
void apply_displacement_bc_csr(int* rowIndex, int* columns, double* values,
    double* F, int n, const DisplacementBC* bc);
void apply_force_bc(double* F, const ForceBC* bc);

#endif
