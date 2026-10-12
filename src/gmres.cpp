#include "krylov_internal.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <new>
#include <stdexcept>

using namespace krylov_detail;
void mvm(Complex*, Complex*);
void mvm_in(Complex*, Complex*);
void mvm_Zneardotx(Complex*, Complex*);

namespace {
Complex phase(Complex value) {
    return value == Complex{} ? Complex{1.0} : value/std::abs(value);
}
void reflect(const Complex* q, Complex* v, int n) {
    const Complex projection = 2.0*dot(q, v, n);
    for (int k = 0; k < n; ++k) v[k] -= q[k]*projection;
}
void reflector(const Complex* v, Complex* q, int first, int n, double norm) {
    std::fill_n(q, n, Complex{});
    // Normalize before adding the phase to avoid v[first] + phase*norm overflow.
    for (int k = first; k < n; ++k) q[k] = v[k]/norm;
    q[first] += phase(v[first]);
    const double qnorm = krylov_norm(q, n);
    for (int k = first; k < n; ++k) q[k] /= qnorm;
}
}

KrylovResult gmres_solve(int n, Complex* b, int maxit, double tolerance,
                        Complex* x, int pc_in, bool verbose) {
    if (pc_in < 0 || pc_in > 2 || !valid_input(n, b, maxit, tolerance, x)) return {};
    if (!finite(b, n) || !finite(x, n)) { KrylovResult r; r.status = KrylovStatus::NonFinite; return r; }
    KrylovResult saved;
    try {
        // Preserve IBC operator selectors: accurate, approximate, near-only.
        const Operator op = pc_in == 0 ? mvm : (pc_in == 1 ? mvm_in : mvm_Zneardotx);
        Context ctx(n, b, x, op, tolerance);
        if (!std::isfinite(ctx.rhs_norm) || !ctx.check(x)) return ctx.finish(KrylovStatus::NonFinite);
        if (ctx.converged()) return ctx.finish(KrylovStatus::Converged);
        if (maxit == 0) return ctx.finish(KrylovStatus::IterationLimit);
        saved = ctx.result;
        const int m = std::min(maxit, n);
        Vector q(storage_size(std::size_t(m)+1, n)), h(storage_size(m, m)), g(std::size_t(m)+1);
        // The extra slot is the zero subdiagonal at the final dimension.
        Vector v(std::size_t(n)+1), product(n), correction(n), candidate(n), base(x, x+n), y;
        std::vector<Rotation> rotations(m);
        reflector(ctx.residual.data(), q.data(), 0, n, ctx.result.residual_norm);
        g[0] = -phase(ctx.residual[0])*ctx.result.residual_norm;
        for (int j = 0; j < m; ++j) {
            const Complex* qj = q.data()+std::size_t(j)*n;
            for (int k = 0; k < n; ++k) v[k] = -2.0*std::conj(qj[j])*qj[k];
            v[j] += 1.0;
            for (int k = j-1; k >= 0; --k) reflect(q.data()+std::size_t(k)*n, v.data(), n);
            const double vnorm = krylov_norm(v.data(), n);
            if (!std::isfinite(vnorm) || vnorm == 0.0) return ctx.finish(KrylovStatus::Breakdown);
            for (int k = 0; k < n; ++k) v[k] /= vnorm;
            if (!ctx.product(v.data(), product.data())) return ctx.finish(KrylovStatus::NonFinite);
            const double before = krylov_norm(product.data(), n);
            std::copy(product.begin(), product.end(), v.begin());
            v[n] = 0.0;
            for (int k = 0; k <= j; ++k) reflect(q.data()+std::size_t(k)*n, v.data(), n);
            const double tail = krylov_norm(v.data()+j+1, n-j-1);
            if (!std::isfinite(before) || !std::isfinite(tail) || !finite(v.data(), n))
                return ctx.finish(KrylovStatus::NonFinite);
            const bool breakdown = j+1 == n || tail <= 16.0*std::numeric_limits<double>::epsilon()*before;
            if (!breakdown) {
                reflector(v.data(), q.data()+std::size_t(j+1)*n, j+1, n, tail);
                v[j+1] = -phase(v[j+1])*tail;
                std::fill(v.begin()+j+2, v.end(), Complex{});
            } else {
                std::fill(v.begin()+j+1, v.end(), Complex{});
            }
            for (int k = 0; k < j; ++k) rotations[k].apply(v[k], v[k+1]);
            if (!rotations[j].make(v[j], v[j+1])) return ctx.finish(KrylovStatus::NonFinite);
            rotations[j].apply(v[j], v[j+1]);
            rotations[j].apply(g[j], g[j+1]);
            for (int k = 0; k <= j; ++k) h[std::size_t(k)*m+j] = v[k];
            if (!finite(v.data(), n+std::size_t(1)) || !finite(g.data(), g.size()))
                return ctx.finish(KrylovStatus::NonFinite);
            ++ctx.result.iterations;
            const double estimate = std::abs(g[j+1])/(ctx.rhs_norm == 0.0 ? 1.0 : ctx.rhs_norm);
            if (verbose) std::printf("GMRES-Iter.: %d, estimated relative residual: %.6e\n", ctx.result.iterations, estimate);
            if (estimate <= tolerance || breakdown || j+1 == m) {
                int usable = j+1;
                while (usable > 0) {
                    const KrylovStatus solve_status = triangular(h, m, g, usable, y);
                    if (solve_status == KrylovStatus::Converged) break;
                    if (solve_status == KrylovStatus::NonFinite) return ctx.finish(solve_status);
                    --usable;
                }
                if (usable == 0) return ctx.finish(KrylovStatus::Breakdown);
                std::fill(correction.begin(), correction.end(), Complex{});
                // Reconstruct the Householder basis combination without storing a second basis.
                for (int k = usable-1; k >= 0; --k) {
                    correction[k] += y[k];
                    reflect(q.data()+std::size_t(k)*n, correction.data(), n);
                }
                for (int k = 0; k < n; ++k) candidate[k] = base[k]+correction[k];
                const double previous_norm = ctx.result.residual_norm;
                const double previous_relative = ctx.result.relative_residual;
                if (!ctx.check(candidate.data())) return ctx.finish(KrylovStatus::NonFinite);
                if (ctx.result.residual_norm <= previous_norm) {
                    std::copy(candidate.begin(), candidate.end(), x);
                    saved = ctx.result;
                } else {
                    ctx.result.residual_norm = previous_norm;
                    ctx.result.relative_residual = previous_relative;
                }
                if (ctx.converged()) return ctx.finish(KrylovStatus::Converged);
                if (breakdown || usable != j+1) return ctx.finish(KrylovStatus::Breakdown);
                // A small recurrence estimate is insufficient: continue if the true residual fails.
            }
        }
        return ctx.finish(KrylovStatus::IterationLimit);
    } catch (const std::bad_alloc&) {
        saved.status = KrylovStatus::AllocationFailure; return saved;
    } catch (const std::length_error&) {
        saved.status = KrylovStatus::AllocationFailure; return saved;
    }
}

int gmres(int n, Complex* b, int maxit, double tolerance, Complex* x, int pc_in) {
    return static_cast<int>(gmres_solve(n, b, maxit, tolerance, x, pc_in, pc_in == 0).status);
}
