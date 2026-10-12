// Small coupled-IBC operator checks. Dense reference storage is test-only.
#include "main.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#ifndef IBC_LEGACY_REFERENCE
#include "mlfma_workspace.h"
#endif
using C=std::complex<double>;
void mvm_Zneardotx(C*,C*);
void mvm_Zfar1121dotx(C*,C*);
void mvm_Zfar1222dotx(C*,C*);
void mvm_Zfar1121dotx_in(C*,C*);
void mvm_Zfar1222dotx_in(C*,C*);
void mvm_in(C*,C*);
void mvm(C*,C*);
void calculate_dir_pw(int*,double,double,double*,int,int);
void calculate_rhsvector(double*);
void mvm_BDPmatrix(C*,C*);
#ifdef IBC_LEGACY_REFERENCE
// The original partition builder writes past its allocation even with one
// thread. A guarded serial reference bypasses it without altering any kernels.
void calculate_loadbalance() {
    delete[] thread_boxes; delete[] thread_boxes_Rad;
    const int first=num_filledboxes[2*num_mlfmalevels-1];
    const int last=num_filledboxes[3*num_mlfmalevels-1]+1;
    thread_boxes=new int[2]{first,last};
    thread_boxes_Rad=new int[2]{first,last};
}
#endif
static void dump(const char* name,const C* data,size_t count) {
    std::ofstream f(name,std::ios::binary);
    f.write(reinterpret_cast<const char*>(data),count*sizeof(C));
    assert(bool(f));
}
static double relative(const std::vector<C>& a,const std::vector<C>& b) {
    assert(a.size()==b.size()); double error=0,norm=0;
    for(size_t i=0;i<a.size();++i) {
        assert(std::isfinite(std::abs(a[i])) && std::isfinite(std::abs(b[i])));
        error+=std::norm(a[i]-b[i]); norm+=std::norm(b[i]);
    }
    return std::sqrt(error/(norm==0 ? 1 : norm));
}
static std::vector<C> explicit_near(const std::vector<C>& x) {
    std::vector<C> y(2*num_rwgs);
    const int first=num_filledboxes[2*num_mlfmalevels-1];
    const int last=num_filledboxes[3*num_mlfmalevels-1]+1;
    for(int box=first;box<last;++box) for(int neighbor=0;neighbor<27;++neighbor) {
        const int basis=boxes_near[27*box+neighbor]; if(basis<0) continue;
        long index=nearindex[27*(box-first)+neighbor];
        for(int r=children[2*box];r<=children[2*box+1];++r)
            for(int c=children[2*basis];c<=children[2*basis+1];++c,++index) {
                const int row=rwg_ind[r],col=rwg_ind[c];
                y[row]+=matrix_Znear11[index]*x[col];
                y[row]+=matrix_Znear12[index]*x[col+num_rwgs];
                y[row+num_rwgs]+=matrix_Znear21[index]*x[col];
                y[row+num_rwgs]+=matrix_Znear22[index]*x[col+num_rwgs];
            }
    }
    return y;
}
#ifndef IBC_LEGACY_REFERENCE
static void check_preconditioner(const std::vector<C>& x) {
    const int saved_pc=val_pc;
    val_pc=0; calculate_blockdiagonalpreconditioner();
    assert(!blockdiagonalmatrix11 && !blockdiagonalmatrix22);
    val_pc=2; calculate_blockdiagonalpreconditioner();
    assert(!blockdiagonalmatrix11 && !blockdiagonalmatrix22);
    val_pc=1; calculate_blockdiagonalpreconditioner();
    std::vector<C> out(x.size()),alias=x,product(x.size());
    mvm_BDPmatrix(const_cast<C*>(x.data()),out.data());
    mvm_BDPmatrix(alias.data(),alias.data());
    assert(relative(alias,out)==0.0);
    const int first=num_filledboxes[2*num_mlfmalevels-1];
    const int last=num_filledboxes[3*num_mlfmalevels-1]+1;
    for(int box=first;box<last;++box) for(int neighbor=0;neighbor<27;++neighbor) {
        if(boxes_near[27*box+neighbor]!=box) continue;
        long index=nearindex[27*(box-first)+neighbor];
        for(int r=children[2*box];r<=children[2*box+1];++r)
            for(int c=children[2*box];c<=children[2*box+1];++c,++index) {
                const int row=rwg_ind[r],col=rwg_ind[c];
                product[row]+=matrix_Znear11[index]*out[col];
                product[row+num_rwgs]+=matrix_Znear22[index]*out[col+num_rwgs];
            }
    }
    const double error=relative(product,x);
    printf("CHECK block_relative_residual %.17g\n",error); assert(error<1e-12);
    val_pc=saved_pc; calculate_blockdiagonalpreconditioner();
    assert(!blockdiagonalmatrix11 && !blockdiagonalmatrix22);
}
static void check_dense(const std::vector<C>& x,const std::vector<C>& accurate) {
    const int first=num_filledboxes[2*num_mlfmalevels-1];
    const int last=num_filledboxes[3*num_mlfmalevels-1]+1;
    assert(last-first<=27 && num_rwgs<=48);
    auto* saved_boxes=boxes_near; auto* saved_index=nearindex;
    auto* m11=matrix_Znear11;auto* m12=matrix_Znear12;
    auto* m21=matrix_Znear21;auto* m22=matrix_Znear22;
    const long saved_size=size_Znear;
    std::vector<int> all_boxes(27*num_filledboxesall,-1);
    std::vector<long> all_index(27*(last-first),-1);
    long offset=0;
    for(int row=first;row<last;++row) for(int col=first;col<last;++col) {
        all_boxes[27*row+col-first]=col;
        all_index[27*(row-first)+col-first]=offset;
        offset+=long(num_children[row])*num_children[col];
    }
    boxes_near=all_boxes.data(); nearindex=all_index.data(); size_Znear=offset;
    calculate_loadbalance(); calculate_Znearmatrix();
    const auto direct=explicit_near(x);
    dump("direct.bin",direct.data(),direct.size());
    const double error=relative(accurate,direct);
    delete[] matrix_Znear11;delete[] matrix_Znear12;
    delete[] matrix_Znear21;delete[] matrix_Znear22;
    boxes_near=saved_boxes;nearindex=saved_index;size_Znear=saved_size;
    matrix_Znear11=m11;matrix_Znear12=m12;matrix_Znear21=m21;matrix_Znear22=m22;
    calculate_loadbalance();
    printf("CHECK mlfma_direct_relative_error %.17g tolerance %.17g\n",error,val_soltol);
    assert(error<=val_soltol);
}
#endif
int main(int argc,char** argv) {
    if(argc!=3) return 2;
    read_simulationinput();
    parallel_omp_on=std::atoi(argv[1]);num_threads=std::atoi(argv[2]);
    generate_edgeconnection();generate_octtree();generate_nearfarboxes();
    generate_Znearcomponent();generate_Zfarcomponents();generate_translationcomponent();
    calculate_loadbalance();calculate_Znearmatrix();calculate_translationoperator();
    aggregation=new C[2*size_agg];disaggregation=new C[2*size_agg];
    const int n=2*num_rwgs;
    std::vector<C> x(n),near(n),far_j(n),far_m(n),approximate(n),accurate(n);
    for(int k=0;k<n;++k) x[k]={std::sin(0.13*k),std::cos(0.17*k)};
    mvm_Zneardotx(x.data(),near.data());
    const double near_error=relative(near,explicit_near(x));
    printf("CHECK near_relative_error %.17g\n",near_error);assert(near_error<1e-12);
    mvm_Zfar1121dotx(x.data(),far_j.data());
    mvm_Zfar1222dotx(x.data(),far_m.data());
    if(val_approx>0) {
        mvm_Zfar1121dotx_in(x.data(),approximate.data());
        mvm_Zfar1222dotx_in(x.data(),approximate.data());
    }
    for(int k=0;k<n;++k) accurate[k]=near[k]+far_j[k]+far_m[k];
    double far_norm=0;for(int k=0;k<n;++k) far_norm+=std::norm(far_j[k])+std::norm(far_m[k]);
    printf("CHECK far_norm %.17g\n",std::sqrt(far_norm));assert(far_norm>0);
    dump("near11.bin",matrix_Znear11,size_Znear);dump("near12.bin",matrix_Znear12,size_Znear);
    dump("near21.bin",matrix_Znear21,size_Znear);dump("near22.bin",matrix_Znear22,size_Znear);
    dump("near.bin",near.data(),n);dump("far_j.bin",far_j.data(),n);
    dump("far_m.bin",far_m.data(),n);dump("approximate.bin",approximate.data(),n);
#ifndef IBC_LEGACY_REFERENCE
    std::vector<C> approx_reference(n),again(n);
    mvm_in(x.data(),approx_reference.data());
    const auto retained=mlfma_workspace().retained_bytes();
    for(int repeat=0;repeat<12;++repeat) {
        std::fill(again.begin(),again.end(),C{});mvm_in(x.data(),again.data());
        assert(relative(again,approx_reference)==0.0);
        assert(mlfma_workspace().retained_bytes()==retained);
    }
    printf("CHECK stable_workspace_bytes %zu\n",retained);
    check_preconditioner(x);check_dense(x,accurate);
    // Independently check full-solver residuals in Python using the materialized
    // tiny accurate operator. This is bounded to 48 unknowns and never production.
    std::vector<C> dense(size_t(n)*n),unit(n),column(n);
    for(int col=0;col<n;++col) {
        std::fill(unit.begin(),unit.end(),C{});unit[col]=1;
        std::fill(column.begin(),column.end(),C{});mvm(unit.data(),column.data());
        for(int row=0;row<n;++row) dense[size_t(row)*n+col]=column[row];
    }
    dump("accurate_matrix.bin",dense.data(),dense.size());
    dir_PW=new double[4];
    calculate_dir_pw(num_source,source_PW[1]-source_PW[0],source_PW[3]-source_PW[2],dir_PW,
                     num_source[0]-1,num_source[1]-1);
    calculate_rhsvector(dir_PW);dump("rhs.bin",vec_RHS,n);
    clear_memory();clear_memory();
    assert(!node_pos && !vec_sol && !matrix_Znear11 && !matrix_Znear12 &&
           !matrix_Znear21 && !matrix_Znear22 && !agg_index && !agg_index_pc &&
           !trans_op && !trans_op_pc && !blockdiagonalmatrix11 && !blockdiagonalmatrix22);
    assert(mlfma_workspace().retained_bytes()==0);
    read_simulationinput();generate_edgeconnection();generate_octtree();clear_memory();clear_memory();
#endif
    puts("PASS coupled IBC operators");
}
