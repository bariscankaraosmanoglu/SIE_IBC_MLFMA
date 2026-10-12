#include "krylov_internal.h"

using namespace krylov_detail;
void mvm_in(Complex*, Complex*);
void pc_solution_in(Complex*, Complex*, int);

KrylovResult fgmres_in_solve(int n, Complex* b, int maxit, int restart,
                          double tolerance, Complex* x, bool verbose) {
    return flexible(n, b, maxit, restart, tolerance, x, mvm_in, pc_solution_in, verbose);
}

int fgmres_in(int n, Complex* b, int maxit, int restart, double tolerance, Complex* x) {
    return static_cast<int>(fgmres_in_solve(n, b, maxit, restart, tolerance, x, false).status);
}
