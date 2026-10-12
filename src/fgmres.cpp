#include "krylov_internal.h"

using namespace krylov_detail;
void mvm(Complex*, Complex*);
void pc_solution(Complex*, Complex*, int);

KrylovResult fgmres_solve(int n, Complex* b, int maxit, int restart,
                          double tolerance, Complex* x, bool verbose) {
    return flexible(n, b, maxit, restart, tolerance, x, mvm, pc_solution, verbose);
}

int fgmres(int n, Complex* b, int maxit, int restart, double tolerance, Complex* x) {
    return static_cast<int>(fgmres_solve(n, b, maxit, restart, tolerance, x, true).status);
}
