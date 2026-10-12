#include "krylov_internal.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <new>
#include <stdexcept>

const char* krylov_status_name(KrylovStatus status) {
    switch (status) {
    case KrylovStatus::Converged: return "converged";
    case KrylovStatus::IterationLimit: return "iteration_limit";
    case KrylovStatus::Breakdown: return "breakdown";
    case KrylovStatus::InvalidInput: return "invalid_input";
    case KrylovStatus::NonFinite: return "nonfinite_arithmetic";
    case KrylovStatus::AllocationFailure: return "allocation_failure";
    case KrylovStatus::ResidualMismatch: return "physical_residual_not_converged";
    }
    return "unknown";
}

double krylov_norm(const std::complex<double>* values, int n) {
    // Scaled sum of squares avoids overflow/underflow from squaring the entries.
    double scale = 0.0, sum = 1.0;
    for (int i = 0; i < n; ++i) {
        for (double v : {values[i].real(), values[i].imag()}) {
            const double a = std::abs(v);
            if (!std::isfinite(a)) return a;
            if (a == 0.0) continue;
            if (scale < a) { const double r = scale / a; sum = 1.0 + sum*r*r; scale = a; }
            else { const double r = a / scale; sum += r*r; }
        }
    }
    return scale * std::sqrt(sum);
}

namespace krylov_detail {
bool finite(const Complex* x, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(x[i].real()) || !std::isfinite(x[i].imag())) return false;
    return true;
}
Complex dot(const Complex* x, const Complex* y, int n) {
    Complex value = 0.0;
    for (int i = 0; i < n; ++i) value += std::conj(x[i])*y[i];
    return value;
}
bool valid_input(int n, const Complex* b, int maxit, double tol, const Complex* x) {
    return n > 0 && b && x && b != x && maxit >= 0 && std::isfinite(tol) && tol >= 0.0;
}
std::size_t storage_size(std::size_t rows, std::size_t columns) {
    if (columns && rows > Vector().max_size()/columns) throw std::length_error("Krylov workspace");
    return rows*columns;
}
bool Rotation::make(Complex a, Complex b) {
    const double rho = std::hypot(std::abs(a), std::abs(b));
    if (!std::isfinite(rho)) return false;
    if (rho == 0.0) { c = 1.0; s = 0.0; }
    else { c = a/rho; s = b/rho; }
    return true;
}
void Rotation::apply(Complex& a, Complex& b) const {
    const Complex old = a;
    a = std::conj(c)*old + std::conj(s)*b;
    b = -s*old + c*b;
}
KrylovStatus triangular(const Vector& h, int stride, const Vector& g, int count, Vector& y) {
    y.assign(g.begin(), g.begin()+count);
    for (int row = count-1; row >= 0; --row) {
        double scale = 0.0;
        for (int col = row; col < count; ++col)
            scale = std::max(scale, std::abs(h[std::size_t(row)*stride+col]));
        const Complex pivot = h[std::size_t(row)*stride+row];
        if (!std::isfinite(scale) || scale == 0.0 ||
            std::abs(pivot) <= std::numeric_limits<double>::epsilon()*scale) return KrylovStatus::Breakdown;
        for (int col = row+1; col < count; ++col)
            y[row] -= h[std::size_t(row)*stride+col]*y[col];
        y[row] /= pivot;
        if (!finite(&y[row], 1)) return KrylovStatus::NonFinite;
    }
    return KrylovStatus::Converged;
}
Context::Context(int dimension, Complex* rhs, Complex* solution, Operator matrix, double tol)
    : n(dimension), b(rhs), x(solution), op(matrix), tolerance(tol),
      rhs_norm(krylov_norm(rhs, dimension)), residual(dimension), work(dimension) {}
bool Context::product(Complex* input, Complex* output) {
    std::fill(output, output+n, Complex{});
    op(input, output);
    ++result.matrix_products;
    return finite(output, n);
}
bool Context::check(Complex* candidate) {
    if (!finite(candidate, n) || !product(candidate, work.data())) return false;
    for (int i = 0; i < n; ++i) work[i] = b[i]-work[i];
    const double norm = krylov_norm(work.data(), n);
    const double relative = norm/(rhs_norm == 0.0 ? 1.0 : rhs_norm);
    if (!std::isfinite(norm) || !std::isfinite(relative)) return false;
    residual = work;
    result.residual_norm = norm;
    result.relative_residual = relative;
    return true;
}
bool Context::converged() const { return result.relative_residual <= tolerance; }
KrylovResult Context::finish(KrylovStatus status) { result.status = status; return result; }

KrylovResult flexible(int n, Complex* b, int maxit, int restart, double tolerance,
                      Complex* x, Operator op, Preconditioner pc, bool verbose) {
    if (!valid_input(n, b, maxit, tolerance, x) || restart <= 0) return {};
    if (!finite(b, n) || !finite(x, n)) { KrylovResult r; r.status = KrylovStatus::NonFinite; return r; }
    KrylovResult saved;
    try {
        Context ctx(n, b, x, op, tolerance);
        if (!std::isfinite(ctx.rhs_norm) || !ctx.check(x)) return ctx.finish(KrylovStatus::NonFinite);
        if (ctx.converged()) return ctx.finish(KrylovStatus::Converged);
        if (maxit == 0) return ctx.finish(KrylovStatus::IterationLimit);
        saved = ctx.result;
        const int m = std::min({n, maxit, restart});
        // Z holds only this restart cycle, independent of the total iteration budget.
        Vector v(storage_size(std::size_t(m)+1, n)), z(storage_size(m, n));
        Vector h(storage_size(std::size_t(m)+1, m)), g(std::size_t(m)+1), y, candidate(n), input(n);
        std::vector<Rotation> rotations(m);
        while (ctx.result.iterations < maxit) {
            std::fill(h.begin(), h.end(), Complex{});
            std::fill(g.begin(), g.end(), Complex{});
            const double beta = ctx.result.residual_norm;
            for (int k = 0; k < n; ++k) v[k] = ctx.residual[k]/beta;
            g[0] = beta;
            int count = 0;
            bool breakdown = false;
            for (int j = 0; j < m && ctx.result.iterations < maxit; ++j) {
                std::copy_n(v.data()+std::size_t(j)*n, n, input.data());
                Complex* zj = z.data()+std::size_t(j)*n;
                std::fill_n(zj, n, Complex{});
                pc(input.data(), zj, n);
                if (!finite(zj, n)) return ctx.finish(KrylovStatus::NonFinite);
                Complex* w = v.data()+std::size_t(j+1)*n;
                if (!ctx.product(zj, w)) return ctx.finish(KrylovStatus::NonFinite);
                const double before = krylov_norm(w, n);
                for (int k = 0; k <= j; ++k) {
                    const Complex* vk = v.data()+std::size_t(k)*n;
                    const Complex proj = dot(vk, w, n);
                    h[std::size_t(k)*m+j] = proj;
                    for (int t = 0; t < n; ++t) w[t] -= proj*vk[t];
                }
                double tail = krylov_norm(w, n);
                // A second pass when cancellation is substantial preserves orthogonality.
                if (tail < 0.5*before) {
                    for (int k = 0; k <= j; ++k) {
                        const Complex* vk = v.data()+std::size_t(k)*n;
                        const Complex proj = dot(vk, w, n);
                        h[std::size_t(k)*m+j] += proj;
                        for (int t = 0; t < n; ++t) w[t] -= proj*vk[t];
                    }
                    tail = krylov_norm(w, n);
                }
                if (!std::isfinite(before) || !std::isfinite(tail))
                    return ctx.finish(KrylovStatus::NonFinite);
                breakdown = tail <= 16.0*std::numeric_limits<double>::epsilon()*before || j+1 == n;
                h[std::size_t(j+1)*m+j] = breakdown ? 0.0 : tail;
                if (!breakdown) for (int k = 0; k < n; ++k) w[k] /= tail;
                for (int k = 0; k < j; ++k)
                    rotations[k].apply(h[std::size_t(k)*m+j], h[std::size_t(k+1)*m+j]);
                if (!rotations[j].make(h[std::size_t(j)*m+j], h[std::size_t(j+1)*m+j]))
                    return ctx.finish(KrylovStatus::NonFinite);
                rotations[j].apply(h[std::size_t(j)*m+j], h[std::size_t(j+1)*m+j]);
                rotations[j].apply(g[j], g[j+1]);
                if (!finite(g.data()+j, 2))
                    return ctx.finish(KrylovStatus::NonFinite);
                for (int k = 0; k <= j+1; ++k)
                    if (!finite(&h[std::size_t(k)*m+j], 1)) return ctx.finish(KrylovStatus::NonFinite);
                ++ctx.result.iterations;
                count = j+1;
                const double estimate = std::abs(g[j+1])/(ctx.rhs_norm == 0.0 ? 1.0 : ctx.rhs_norm);
                if (verbose) std::printf("FGMRES-Iter.: %d, estimated relative residual: %.6e\n", ctx.result.iterations, estimate);
                if (estimate <= tolerance || breakdown) break;
            }
            int usable = count;
            while (usable > 0) {
                const KrylovStatus solve_status = triangular(h, m, g, usable, y);
                if (solve_status == KrylovStatus::Converged) break;
                if (solve_status == KrylovStatus::NonFinite) return ctx.finish(solve_status);
                --usable;
            }
            if (usable == 0) return ctx.finish(KrylovStatus::Breakdown);
            breakdown = breakdown || usable != count;
            std::copy_n(x, n, candidate.data());
            for (int j = 0; j < usable; ++j)
                for (int k = 0; k < n; ++k) candidate[k] += y[j]*z[std::size_t(j)*n+k];
            if (!ctx.check(candidate.data())) return ctx.finish(KrylovStatus::NonFinite);
            std::copy(candidate.begin(), candidate.end(), x);
            saved = ctx.result;
            if (ctx.converged()) return ctx.finish(KrylovStatus::Converged);
            if (breakdown) return ctx.finish(KrylovStatus::Breakdown);
            // Restart from the explicitly computed residual, even if the estimate was small.
        }
        return ctx.finish(KrylovStatus::IterationLimit);
    } catch (const std::bad_alloc&) {
        saved.status = KrylovStatus::AllocationFailure; return saved;
    } catch (const std::length_error&) {
        saved.status = KrylovStatus::AllocationFailure; return saved;
    }
}
}
