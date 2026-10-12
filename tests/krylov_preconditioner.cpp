#include "krylov.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
using C = std::complex<double>;
int val_pc = 2, val_itermax_pc = 1;
static bool bad_operator;
void mvm(C* x, C* y) { y[0]=C{0,1}*x[0]+x[1]; y[1]=x[0]+C{0,1}*x[1]; }
void mvm_in(C* x, C* y) {
    mvm(x,y);
    if (bad_operator) y[0] = std::numeric_limits<double>::infinity();
}
static int near_products;
void mvm_Zneardotx(C* x,C* y) { ++near_products; mvm(x,y); }

void pc_solution(C*, C*, int);
int main() {
    C b[2]={1.0,0.0}, x[2]{};
    pc_solution(b,x,2);
    assert(near_products>0);
    assert(std::isfinite(std::abs(x[0])) && std::isfinite(std::abs(x[1])));
    assert(std::abs(x[0])+std::abs(x[1]) > 0);
    x[0]=x[1]=0.0;
    const auto result=fgmres_solve(2,b,2,2,1e-12,x);
    assert(result.status==KrylovStatus::Converged && result.relative_residual<=1e-12);
    std::printf("CHECK finite_limited_inner_solve outer_residual=%.17g\n",result.relative_residual);
    val_itermax_pc=-1;
    x[0]=x[1]=0.0;
    const auto invalid=fgmres_solve(2,b,2,2,1e-12,x);
    assert(invalid.status==KrylovStatus::NonFinite && x[0]==C{} && x[1]==C{});
    val_itermax_pc=2; bad_operator=true;
    const auto nonfinite=fgmres_solve(2,b,2,2,1e-12,x);
    assert(nonfinite.status==KrylovStatus::NonFinite && x[0]==C{} && x[1]==C{});
    std::puts("PASS inner preconditioner status checks");
}
