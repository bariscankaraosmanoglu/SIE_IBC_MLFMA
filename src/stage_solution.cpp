#include "solver_globals.h"
#include "krylov.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <complex>
#include <chrono>
#include <math.h>
#include "matrixBlock.h"
#include <iostream>

using namespace std;
using namespace std::chrono;


extern int num_rwgs;
extern int val_itermax;
extern long int size_agg;
extern complex<double> *aggregation;
extern complex<double> *disaggregation;
extern complex<double> *vec_RHS;
extern complex<double> *vec_sol;
extern double val_soltol;
extern double *dir_PW;
extern int *num_source;
extern int val_pol;
extern int val_pc;
extern double dir_ff;

extern int *num_filledboxes;
extern int num_mlfmalevels;


void calculate_dir_pw(int num_source[], double delta_theta, double delta_phi, double dir_PW[], int ind_theta, int ind_phi);
void calculate_initguess(int ind_sol, complex <double> val_RHS_ratio, complex <double> vec_RHS_prev[]);
void calculate_rhsvector(double dir_PW[]);
void mvm_physical(complex<double> x[], complex<double> y[]);
void copy_rhsvector(complex <double> vec_RHS_prev[]);
void calculate_farfield(int ind_sol);

void mvm_BDPmatrix(std::complex<double> x[], std::complex<double> y[]);


int stage_solution(){

    const int num_unknowns=2*num_rwgs;
    int itermax=val_itermax;
    int iterrestart=val_itermax;
    double soltolerance=val_soltol;
    
    double delta_theta=(source_PW[1]-source_PW[0]);
    double delta_phi=(source_PW[3]-source_PW[2]);
    vector<complex<double>> vec_RHS_prev(num_unknowns), physical_rhs(num_unknowns), physical_product(num_unknowns);
    ofstream status_file("result_solver_status.txt");
    if (!status_file) fprintf(stderr, "Cannot open result_solver_status.txt\n");
    status_file << setprecision(17)
                << "# illumination status krylov_status iterations matrix_products krylov_true_relative_residual physical_true_relative_residual tolerance\n";
    int exit_status = status_file ? 0 : 1;
    dir_PW = new double[4];
    aggregation=new complex<double>[size_agg*2];
    disaggregation=new complex<double>[size_agg*2];
    vec_sol=new complex<double>[num_unknowns];
    for(int ind0=0;ind0<num_unknowns;ind0++){
        vec_sol[ind0]=complex<double>(0.0, 0.0);
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    int ind_sol=0;
    for(int ind_theta=0; ind_theta<num_source[0]; ind_theta++){
        for(int ind_phi=0; ind_phi<num_source[1]; ind_phi++){
        calculate_dir_pw(num_source, delta_theta, delta_phi, dir_PW, ind_theta, ind_phi);
        calculate_rhsvector(dir_PW);
        copy_n(vec_RHS, num_unknowns, physical_rhs.data());
        if (val_pc==1){mvm_BDPmatrix(vec_RHS, vec_RHS);}
        calculate_initguess(ind_sol, complex<double>{}, vec_RHS_prev.data());

        KrylovResult result;
        if((val_pc==0)||(val_pc==1)) result = gmres_solve(num_unknowns,vec_RHS,itermax,soltolerance,vec_sol,0,true);
        if(val_pc==2) result = fgmres_solve(num_unknowns,vec_RHS,itermax,iterrestart,soltolerance,vec_sol,true);

        // Check b - A*x using the original RHS and accurate, unpreconditioned operator.
        fill(physical_product.begin(), physical_product.end(), complex<double>{});
        mvm_physical(vec_sol, physical_product.data());
        for (int k = 0; k < num_unknowns; ++k) physical_product[k] = physical_rhs[k]-physical_product[k];
        const double rhs_norm = krylov_norm(physical_rhs.data(), num_unknowns);
        const double physical_relative = krylov_norm(physical_product.data(), num_unknowns)/(rhs_norm == 0.0 ? 1.0 : rhs_norm);
        KrylovStatus status = result.status;
        if (!isfinite(physical_relative)) status = KrylovStatus::NonFinite;
        else if (status == KrylovStatus::Converged && physical_relative > soltolerance)
            status = KrylovStatus::ResidualMismatch;
        if (status != KrylovStatus::Converged) exit_status = 1;
        printf("Solver illumination %d: %s, iterations=%d, true relative residual=%.6e, physical relative residual=%.6e, tolerance=%.6e\n",
               ind_sol, krylov_status_name(status), result.iterations, result.relative_residual, physical_relative, soltolerance);
        status_file << ind_sol << ' ' << krylov_status_name(status) << ' ' << krylov_status_name(result.status)
                    << ' ' << result.iterations << ' ' << result.matrix_products << ' ' << result.relative_residual
                    << ' ' << physical_relative << ' ' << soltolerance << '\n';
        status_file.flush();

        copy_rhsvector(vec_RHS_prev.data());
        calculate_farfield(ind_sol);

        delete[] vec_RHS;
        vec_RHS = nullptr;
        ind_sol++;
        }
    }



    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    printf("Solution Part Takes %f s \n", time_span.count());

    std::cout.flush();
    return exit_status;
}


//
// Plane wave illumination and polarization directions are calculated from input PW angles
void calculate_dir_pw(int num_source[], double delta_theta, double delta_phi, double dir_PW[], int ind_theta, int ind_phi){
    if (num_source[0]>1){
        delta_theta/=(num_source[0]-1);
    }
    if (num_source[1]>1){
        delta_phi/=(num_source[1]-1);
    }

    dir_PW[0]=source_PW[0]+delta_theta*ind_theta;
    dir_PW[1]=source_PW[2]+delta_phi*ind_phi;
    if (val_pol==0){
        dir_PW[2]=fmod(dir_PW[0]+90, 180);
        dir_PW[3]=dir_PW[1];
    } else{
        dir_PW[2]=90.0;
        dir_PW[3]=fmod(dir_PW[1]+90, 360);
    }

}

//
// Initial guess vector for the next solution is calculated via present and previous rhs vectors
void calculate_initguess(int ind_sol, complex <double> val_RHS_ratio, complex <double> vec_RHS_prev[]){
    if (ind_sol>0){
        for(int ind_rwg=0;ind_rwg<(2*num_rwgs);ind_rwg++){
            if (vec_RHS_prev[ind_rwg] == complex<double>{}) continue;
            val_RHS_ratio=vec_RHS[ind_rwg]/vec_RHS_prev[ind_rwg];
            if (abs(val_RHS_ratio)<2){
                vec_sol[ind_rwg]*=val_RHS_ratio;
            }
        }
    }
}

//
// Copy the present rhs vector for future calculations
void copy_rhsvector(complex <double> vec_RHS_prev[]){
    for(int ind_rwg=0;ind_rwg<(2*num_rwgs);ind_rwg++){
        vec_RHS_prev[ind_rwg]=vec_RHS[ind_rwg];
    }
}