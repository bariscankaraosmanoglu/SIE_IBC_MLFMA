#ifndef KRYLOV_H
#define KRYLOV_H

#include <complex>
#include <limits>

enum class KrylovStatus {
    Converged = 0, IterationLimit, Breakdown, InvalidInput, NonFinite,
    AllocationFailure, ResidualMismatch
};

struct KrylovResult {
    KrylovStatus status = KrylovStatus::InvalidInput;
    int iterations = 0;
    int matrix_products = 0;
    double residual_norm = std::numeric_limits<double>::infinity();
    double relative_residual = std::numeric_limits<double>::infinity();
};

const char* krylov_status_name(KrylovStatus status);
double krylov_norm(const std::complex<double>* values, int n);

// Relative to ||b||; for b == 0 the tolerance and reported residual are absolute.
// On arithmetic failure x retains the last iterate with a checked finite residual.
KrylovResult gmres_solve(int n, std::complex<double>* b, int itermax,
                        double tolerance, std::complex<double>* x, int pc_in,
                        bool verbose = false);
KrylovResult fgmres_solve(int n, std::complex<double>* b, int itermax, int restart,
                         double tolerance, std::complex<double>* x, bool verbose = false);
KrylovResult fgmres_in_solve(int n, std::complex<double>* b, int itermax, int restart,
                            double tolerance, std::complex<double>* x, bool verbose = false);

// Compatibility entry points: zero means verified convergence, otherwise status.
int gmres(int n, std::complex<double>* b, int itermax, double tolerance,
          std::complex<double>* x, int pc_in);
int fgmres(int n, std::complex<double>* b, int itermax, int restart,
           double tolerance, std::complex<double>* x);
int fgmres_in(int n, std::complex<double>* b, int itermax, int restart,
              double tolerance, std::complex<double>* x);

#endif
