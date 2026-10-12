// Dense operators and known solutions independent of the MoM/MLFMA kernels.
#include "krylov_internal.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdio>
#include <string>

using C = std::complex<double>;
using V = std::vector<C>;
static int dimension, pc_calls, products;
static V matrix;
static bool varying_pc, bad_pc, bad_operator;
int parallel_omp_on = 0;
int val_pc = 0;

void mvm(C* x, C* y) {
    ++products;
    for (int row = 0; row < dimension; ++row) {
        y[row] = 0.0;
        for (int col = 0; col < dimension; ++col) y[row] += matrix[row*dimension+col]*x[col];
    }
    if (bad_operator) y[0] = std::numeric_limits<double>::quiet_NaN();
}
void mvm_in(C* x, C* y) { mvm(x, y); }
void mvm_Zneardotx(C* x,C* y) { mvm(x,y); }

void pc_solution(C* b, C* x, int n) {
    ++pc_calls;
    const C factor = varying_pc ? C{0.7+0.15*(pc_calls%3), 0.05*(pc_calls%2)} : C{1.0};
    for (int k = 0; k < n; ++k) x[k] = factor*b[k];
    if (bad_pc) x[0] = std::numeric_limits<double>::quiet_NaN();
}
void pc_solution_in(C* b, C* x, int n) { pc_solution(b, x, n); }

static double norm(const V& values) {
    long double result = 0;
    for (C value : values) result = std::hypot(result, std::hypot((long double)value.real(), (long double)value.imag()));
    return static_cast<double>(result);
}
static double residual(const V& b, const V& x) {
    V r(b);
    for (int row = 0; row < dimension; ++row) {
        C sum{};
        for (int col = 0; col < dimension; ++col) sum += matrix[row*dimension+col]*x[col];
        r[row] -= sum;
    }
    const double nb = norm(b);
    return norm(r)/(nb == 0.0 ? 1.0 : nb);
}
static void setup(int n, V a) {
    dimension = n; matrix = std::move(a); products = pc_calls = 0;
    varying_pc = bad_pc = bad_operator = false;
}
static KrylovResult solve(int method, V& b, V& x, int maxit, int restart = INT_MAX, double tol = 1e-12) {
    if (method == 0) return gmres_solve(dimension, b.data(), maxit, tol, x.data(), 0);
    if (method == 1) return fgmres_solve(dimension, b.data(), maxit, restart, tol, x.data());
    return fgmres_in_solve(dimension, b.data(), maxit, restart, tol, x.data());
}
static void check(const char* name, const KrylovResult& result, const V& b, const V& x, KrylovStatus status) {
    assert(result.status == status);
    assert(krylov_detail::finite(x.data(), x.size()));
    const double actual = residual(b, x);
    assert(std::isfinite(result.relative_residual));
    assert(std::abs(result.relative_residual-actual) <= 2e-14*std::max(1.0, actual));
    assert(result.matrix_products == products);
    if (status == KrylovStatus::Converged) assert(actual <= 1e-12);
    std::printf("CHECK %s status=%s iterations=%d actual_residual=%.17g\n",
                name, krylov_status_name(result.status), result.iterations, actual);
}
static void rotations() {
    for (double scale : {1e-300, 1.0, 1e300}) {
        for (auto pair : {std::pair<C,C>{{0,scale},{scale,0}}, {{scale,2*scale},{-2*scale,scale}},
                          {{0,0},{scale,-scale}}, {{scale,scale},{0,0}}, {{0,0},{0,0}}}) {
            auto a = pair.first, b = pair.second;
            krylov_detail::Rotation rotation;
            assert(rotation.make(a,b));
            rotation.apply(a,b);
            assert(std::abs(b)/scale <= 2e-15);
            assert(std::abs(std::imag(a))/scale <= 2e-15);
            assert(std::abs(std::norm(rotation.c)+std::norm(rotation.s)-1) <= 1e-15);
        }
    }
    assert(krylov_norm(V{C{1e300,0}, C{0,1e300}}.data(), 2) > 1e300);
    assert(krylov_norm(V{C{1e-300,0}, C{0,1e-300}}.data(), 2) > 1e-300);
    bool overflow_checked = false;
    try { krylov_detail::storage_size(std::numeric_limits<std::size_t>::max(), 2); }
    catch (const std::length_error&) { overflow_checked = true; }
    assert(overflow_checked);
}

int main() {
    rotations();
    for (int method = 0; method < 3; ++method) {
        setup(1, {{2,-3}}); V b{{1,4}}, x(1);
        check("one_unknown", solve(method,b,x,1), b,x,KrylovStatus::Converged);
        setup(2, {{0,1},{1,0},{1,0},{0,1}}); b={1.0,0.0}; x.assign(2,{});
        check("complex_cancellation", solve(method,b,x,2), b,x,KrylovStatus::Converged);
        assert(std::abs(x[0]-C{0,-0.5}) < 1e-12 && std::abs(x[1]-C{0.5,0}) < 1e-12);
        x.assign(2,{}); products=0;
        check("iteration_limit", solve(method,b,x,1), b,x,KrylovStatus::IterationLimit);
        assert(std::abs(residual(b,x)-std::sqrt(0.5)) < 1e-14);
        x.assign(2,{}); products=0;
        check("zero_budget", solve(method,b,x,0), b,x,KrylovStatus::IterationLimit);
        x={C{0,-0.5},C{0.5,0}}; products=0;
        const auto exact = solve(method,b,x,0);
        check("exact_initial_guess", exact,b,x,KrylovStatus::Converged);
        assert(exact.iterations == 0);
        b.assign(2,{}); x.assign(2,{}); products=0;
        check("zero_rhs", solve(method,b,x,2), b,x,KrylovStatus::Converged);
        x={C{1,2},C{3,4}}; products=0;
        check("zero_rhs_nonzero_guess", solve(method,b,x,2), b,x,KrylovStatus::Converged);
        setup(2, {1.0,0.0,0.0,1.0}); b={{1,2},{3,-4}}; x.assign(2,{});
        const auto happy = solve(method,b,x,99);
        check("identity_happy_breakdown", happy,b,x,KrylovStatus::Converged);
        assert(happy.iterations == 1);
        setup(2, V(4)); b={1.0,2.0}; x.assign(2,{});
        check("zero_operator", solve(method,b,x,2), b,x,KrylovStatus::Breakdown);
        setup(2, {1.0,0.0,0.0,0.0}); b={1.0,0.0}; x.assign(2,{});
        check("singular_consistent", solve(method,b,x,2), b,x,KrylovStatus::Converged);
        b={1.0,1.0}; x.assign(2,{}); products=0;
        check("singular_inconsistent", solve(method,b,x,2), b,x,KrylovStatus::Breakdown);
        for (double scale : {1e-300,1e-200,1e200,1e300}) {
            setup(2, {{0,scale},{scale,0},{scale,0},{0,scale}});
            b={scale,0.0}; x.assign(2,{});
            check("scaled_complex_system",solve(method,b,x,2),b,x,KrylovStatus::Converged);
        }
        // Known solution for a nonnormal complex matrix. Full-space completion
        // exercises the n boundary; restarting exercises cycle accounting.
        const int n = 12; V a(n*n), known(n);
        for (int row = 0; row < n; ++row) {
            known[row] = C{std::sin(row+0.3), std::cos(row+0.2)};
            a[row*n+row] = C{2+0.05*row, 0.4};
            if (row+1<n) a[row*n+row+1] = C{0.7,-0.3};
        }
        setup(n,a); b.assign(n,{}); mvm(known.data(),b.data()); products=0; x.assign(n,{});
        varying_pc = method != 0;
        check("nonnormal_complex",solve(method,b,x,n),b,x,KrylovStatus::Converged);
        V diff(x); for (int k=0;k<n;++k) diff[k] -= known[k];
        assert(norm(diff)/norm(known) < 1e-11);
        if (method != 0) {
            x.assign(n,{}); products=pc_calls=0;
            const auto restarted = solve(method,b,x,100,2);
            check("flexible_restarts",restarted,b,x,KrylovStatus::Converged);
            assert(restarted.iterations > 2 && restarted.iterations <= 100);
            x.assign(n,{}); products=pc_calls=0;
            check("restart_budget_boundary",solve(method,b,x,3,2),b,x,KrylovStatus::IterationLimit);
            assert(solve(method,b,x,1,0).status == KrylovStatus::InvalidInput);
            assert(solve(method,b,x,1,-1).status == KrylovStatus::InvalidInput);
        }
        assert(solve(method,b,x,-1).status == KrylovStatus::InvalidInput);
        assert(solve(method,b,x,1,1,-0.1).status == KrylovStatus::InvalidInput);
        assert(solve(method,b,x,1,1,std::numeric_limits<double>::quiet_NaN()).status == KrylovStatus::InvalidInput);
        dimension=0;
        assert(solve(method,b,x,1).status == KrylovStatus::InvalidInput); dimension=n;
        x.assign(n,{}); bad_operator=true;
        assert(solve(method,b,x,1).status == KrylovStatus::NonFinite); bad_operator=false;
        b[0]=std::numeric_limits<double>::infinity();
        assert(solve(method,b,x,1).status == KrylovStatus::NonFinite);
        b[0]=1.0; x[0]=std::numeric_limits<double>::quiet_NaN();
        assert(solve(method,b,x,1).status == KrylovStatus::NonFinite);
        x.assign(n,{});
        if (method != 0) {
            bad_pc=true;
            const auto result = solve(method,b,x,1);
            assert(result.status == KrylovStatus::NonFinite && norm(x) == 0.0);
            bad_pc=false;
        }
        // A recurrence reaching zero cannot override the explicitly checked residual.
        setup(n,a); b.assign(n,{}); mvm(known.data(),b.data()); products=0; x.assign(n,{});
        const auto strict = solve(method,b,x,n,INT_MAX,1e-30);
        assert(strict.status != KrylovStatus::Converged && strict.relative_residual > 1e-30);
        check("unattainable_tolerance",strict,b,x,KrylovStatus::Breakdown);
        // Repeated initial convergence, invalid input, and breakdown paths expose leaks.
        for (int repeat=0;repeat<100;++repeat) {
            setup(2,{1.0,0.0,0.0,1.0}); b.assign(2,{}); x.assign(2,{});
            assert(solve(method,b,x,2).status == KrylovStatus::Converged);
            assert(solve(method,b,x,-1).status == KrylovStatus::InvalidInput);
            matrix.assign(4,{}); b[0]=1.0;
            assert(solve(method,b,x,2).status == KrylovStatus::Breakdown);
        }
    }
    std::puts("PASS all Krylov regression checks");
}
