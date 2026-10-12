#ifndef KRYLOV_INTERNAL_H
#define KRYLOV_INTERNAL_H

#include "krylov.h"
#include <vector>

namespace krylov_detail {
using Complex = std::complex<double>;
using Vector = std::vector<Complex>;
using Operator = void (*)(Complex*, Complex*);
using Preconditioner = void (*)(Complex*, Complex*, int);

bool finite(const Complex* x, std::size_t n);
Complex dot(const Complex* x, const Complex* y, int n);
bool valid_input(int n, const Complex* b, int maxit, double tol, const Complex* x);
std::size_t storage_size(std::size_t rows, std::size_t columns);

// U = [conj(c), conj(s); -s, c], U*[a,b] = [rho,0].
struct Rotation {
    Complex c = 1.0, s = 0.0;
    bool make(Complex a, Complex b);
    void apply(Complex& a, Complex& b) const;
};

KrylovStatus triangular(const Vector& h, int stride, const Vector& g, int count, Vector& y);

struct Context {
    int n;
    Complex *b, *x;
    Operator op;
    double tolerance, rhs_norm;
    KrylovResult result;
    Vector residual, work;
    Context(int n, Complex* b, Complex* x, Operator op, double tolerance);
    bool product(Complex* input, Complex* output);
    bool check(Complex* candidate);
    bool converged() const;
    KrylovResult finish(KrylovStatus status);
};

KrylovResult flexible(int n, Complex* b, int maxit, int restart, double tolerance,
                      Complex* x, Operator op, Preconditioner pc, bool verbose);
}
#endif
