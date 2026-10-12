#include "krylov.h"
#include <algorithm>
#include <cstdio>
#include <limits>

using namespace std;
extern int val_pc;
extern int val_itermax_pc;

void pc_solution(complex<double>* b, complex<double>* x, int n) {
    if (val_pc != 2) {
        copy_n(b, n, x);
        return;
    }
    fill_n(x, n, complex<double>{});
    // Preserve the IBC nesting and its original inner tolerance: approximate
    // FGMRES, itself preconditioned by near-only GMRES in pc_solution_in().
    const KrylovResult result = fgmres_in_solve(n,b,val_itermax_pc,val_itermax_pc,0.5,x);
    if (result.status == KrylovStatus::Converged || result.status == KrylovStatus::IterationLimit) return;
    fprintf(stderr, "Inner FGMRES failed: %s, iterations=%d, true relative residual=%.6e\n",
            krylov_status_name(result.status), result.iterations, result.relative_residual);
    // The callback cannot return a status. The outer solver detects this marker
    // before applying its operator and retains its last checked finite iterate.
    fill_n(x, n, complex<double>{numeric_limits<double>::quiet_NaN(), 0.0});
}
