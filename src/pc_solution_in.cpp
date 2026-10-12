#include "krylov.h"
#include <algorithm>
#include <cstdio>
#include <limits>

using namespace std;

void pc_solution_in(complex<double>* b, complex<double>* x, int n) {
    // Both current blocks need initialization. Keep the IBC near-only operator,
    // original 20-iteration budget and 0.4 tolerance of this innermost solve.
    fill_n(x, n, complex<double>{});
    const KrylovResult result = gmres_solve(n,b,20,0.4,x,2);
    if (result.status == KrylovStatus::Converged || result.status == KrylovStatus::IterationLimit) return;
    fprintf(stderr, "Innermost GMRES failed: %s, iterations=%d, true relative residual=%.6e\n",
            krylov_status_name(result.status), result.iterations, result.relative_residual);
    fill_n(x, n, complex<double>{numeric_limits<double>::quiet_NaN(), 0.0});
}
