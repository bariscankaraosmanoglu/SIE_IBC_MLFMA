// Manufactured solutions exercise complex solves independently of factorization.
#include "matrixBlock.h"
#ifndef MLFMA_PIVOTED_BLOCK_LU
#include "matrix_inverse.h"
#endif
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using C = std::complex<double>;
using Values = std::vector<C>;
static bool finite(const C& x) { return std::isfinite(x.real()) && std::isfinite(x.imag()); }
static Values product(const Values& a, const Values& x) {
    Values b(x.size());
    for (std::size_t i=0;i<x.size();++i)
        for (std::size_t j=0;j<x.size();++j) b[i]+=a[i*x.size()+j]*x[j];
    return b;
}
static double error(const Values& a, const Values& b) {
    double scale=0,difference=0;
    for (std::size_t i=0;i<a.size();++i) {
        if (!finite(b[i])) return INFINITY;
        scale=std::max(scale,std::abs(a[i]));difference=std::max(difference,std::abs(a[i]-b[i]));
    }
    return difference/(scale ? scale : 1.0);
}
static Values solve(const Values& a, const Values& b, bool repeat) {
    matrixBlock block;
    const int n=static_cast<int>(b.size());
    block.initmatrixBlock(n,n);
#ifdef MLFMA_PIVOTED_BLOCK_LU
    if (repeat) {
        block.assignvalue(0,0,C(2,1));
        block.initmatrixBlock(1,1);block.assignvalue(0,0,C(2,1));block.factorize();
        C rhs(2,1);block.solve_inplace(&rhs);
        if (std::abs(rhs-C(1,0))>1e-14) throw std::runtime_error("reinitialization solve failed");
        block.initmatrixBlock(n,n);
    }
#else
    (void)repeat;
#endif
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) block.assignvalue(i,j,a[i*n+j]);
#ifdef MLFMA_PIVOTED_BLOCK_LU
    block.factorize();
    Values x=b;block.solve_inplace(x.data());return x;
#else
    auto inv=matrix_inverse(block.submatrix,n);
    Values x(n);
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) x[i]+=inv[i][j]*b[j];
    for (int i=0;i<n;++i) delete[] inv[i];
    delete[] inv;return x;
#endif
}
int main(int argc,char** argv) {
    if (argc!=2) return 2;
    const std::string mode=argv[1];
    if (mode=="empty") {
#ifdef MLFMA_PIVOTED_BLOCK_LU
        matrixBlock empty;C rhs(1,0);
        bool rejected=false;
        try { empty.solve_inplace(&rhs); } catch (const std::exception&) { rejected=true; }
        const bool good=rejected && !empty.ready() && empty.factors.empty();
        std::printf("CHECK empty_storage %d\n",good);return good ? 0 : 1;
#else
        std::puts("CHECK baseline empty block has no safe initialized state");return 1;
#endif
    }
    Values a;int n=3;
    if (mode=="later-pivot" || mode=="repeated") a={1,1,0,1,1,1,0,1,1};
    else if (mode=="first-pivot") { n=2;a={0,C(1,1),C(2,-1),3}; }
    else if (mode=="singular") { n=2;a={1,1,1,1}; }
    else if (mode=="near-singular") { n=2;a={1,1,1,1+1e-10}; }
    else if (mode=="nonfinite") { n=2;a={1,2,3,C(INFINITY,0)}; }
    else {
        n=32;a.resize(n*n);
        for (int i=0;i<n;++i) for(int j=0;j<n;++j)
            a[i*n+j]=i==j ? C(40+i,3-i*0.1) : C(0.2*std::sin(i+2*j),0.15*std::cos(3*i-j));
        if (mode=="scaled-small") for (auto& x:a) x*=1e-150;
        else if (mode=="scaled-large") for (auto& x:a) x*=1e150;
        else if (mode!="complex") return 2;
    }
    Values expected(n);
    for(int i=0;i<n;++i) expected[i]=C(0.5+0.1*i,-0.2+0.03*i);
    const Values rhs=product(a,expected);
    const bool should_reject=mode=="singular" || mode=="nonfinite";
    try {
        const Values x=solve(a,rhs,mode=="repeated");
        if (should_reject) { std::puts("CHECK block was not rejected");return 1; }
        const double forward=error(expected,x),backward=error(rhs,product(a,x));
        std::printf("CHECK forward_error %.17g backward_error %.17g\n",forward,backward);
        return backward<=1e-12 && (mode=="near-singular" || forward<=1e-12) ? 0 : 1;
    } catch(const std::exception& e) {
        std::printf("CHECK rejection %s\n",e.what());return should_reject ? 0 : 1;
    }
}
