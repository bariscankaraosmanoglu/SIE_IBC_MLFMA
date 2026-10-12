// Exercise the actual stage/status plumbing with a known two-unknown operator.
#include "main.h"
#include "krylov.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <string>

static int fields;
static const complex<double> a[4]={{0,1},1.0,1.0,{0,1}};
void calculate_rhsvector(double*) { vec_RHS=new complex<double>[2]{1.0,1.0}; }
void mvm_physical(complex<double>* x,complex<double>* y) {
    y[0]=a[0]*x[0]+a[1]*x[1]; y[1]=a[2]*x[0]+a[3]*x[1];
}
void mvm_BDPmatrix(complex<double>* x,complex<double>* y) {
    y[0]=1e6*x[0]; y[1]=x[1];
}
void mvm(complex<double>* x,complex<double>* y) {
    mvm_physical(x,y); if(val_pc==1) mvm_BDPmatrix(y,y);
}
void mvm_in(complex<double>* x,complex<double>* y) { mvm(x,y); }
void mvm_Zneardotx(complex<double>* x,complex<double>* y) { mvm(x,y); }

void pc_solution(complex<double>* b,complex<double>* x,int n) { std::copy_n(b,n,x); }
void calculate_farfield(int) {
    ++fields;
    assert(isfinite(abs(vec_sol[0])) && isfinite(abs(vec_sol[1])));
}
void calculate_farfield_mlfma(int ind,complex<double>*) { calculate_farfield(ind); }

static void run(const char* name,int pc,int maxit,int illuminations,int expected_exit,const char* expected_status) {
    num_rwgs=1; val_pc=pc; val_itermax=maxit; val_soltol=1e-3;
    size_agg=0; num_source[0]=1; num_source[1]=illuminations;
    source_PW[0]=source_PW[1]=37; source_PW[2]=source_PW[3]=19;
    fields=0;
    assert(stage_solution()==expected_exit);
    assert(fields==illuminations && vec_RHS==nullptr);
    ifstream status_file("result_solver_status.txt");
    string line; getline(status_file,line);
    for(int k=0;k<illuminations;++k) {
        int ind,iterations,products; string status,krylov_status; double residual,physical,tolerance;
        assert(bool(status_file>>ind>>status>>krylov_status>>iterations>>products>>residual>>physical>>tolerance));
        assert(ind==k && status==expected_status && isfinite(residual) && isfinite(physical));
        complex<double> output[2]; mvm_physical(vec_sol,output);
        double independently_checked=hypot(abs(1.0-output[0]),abs(1.0-output[1]))/sqrt(2.0);
        assert(abs(physical-independently_checked)<1e-12);
        if (status=="converged") assert(physical<=tolerance);
        if (status=="physical_residual_not_converged") {
            assert(krylov_status=="converged" && residual<=tolerance && physical>tolerance);
        }
        if (k>0 && status=="converged") assert(iterations==0);
    }
    assert(!(status_file>>line));
    status_file.close();
    assert(std::rename("result_solver_status.txt",(string(name)+"-status.txt").c_str())==0);
    delete[] dir_PW; delete[] aggregation; delete[] disaggregation; delete[] vec_sol;
    dir_PW=nullptr; aggregation=disaggregation=vec_sol=nullptr;
}
int main() {
    run("converged",0,2,1,0,"converged");
    // A strongly weighted preconditioned residual is small after one step,
    // while the physical equation is still wrong by order one.
    run("physical-mismatch",1,2,1,1,"physical_residual_not_converged");
    run("limit",0,0,1,1,"iteration_limit");
    run("multiple",0,2,2,0,"converged");
    // num_source borrows fixed storage in main.h.
    puts("PASS stage physical residual and termination checks");
}
