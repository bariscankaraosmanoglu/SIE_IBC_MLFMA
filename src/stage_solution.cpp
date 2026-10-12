#include <stdio.h>
#include <complex>
#include <chrono>
#include <math.h>
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
extern double source_PW[];
extern double *dir_PW;
extern int *num_source;
extern int val_pol;
extern int val_pc;

void calculate_dir_pw(int num_source[], double delta_theta, double delta_phi, double dir_PW[], int ind_theta, int ind_phi);
void calculate_initguess(int ind_sol, complex <double> val_RHS_ratio, complex <double> vec_RHS_prev[]);
void calculate_rhsvector(double dir_PW[]);
void mvm(complex<double> x[], complex<double> y[]);
int fgmres(int size_vecRHS,complex<double> vec_RHS[],int itermax,int iterrestart,double soltolerance,complex<double> vec_sol[]);
int gmres(int size_vecRHS,complex<double> vec_RHS[],int itermax,double soltolerance,complex<double> vec_sol[], int pc_in);
void copy_rhsvector(complex <double> vec_RHS_prev[]);
void calculate_farfield(int ind_sol);

void mvm(complex<double> x[], complex<double> y[]);
void mvm_Zneardotx(complex<double> x[], complex<double> y[]);
void mvm_Zfar1121dotx(complex<double> x[], complex<double> y[]);
void mvm_Zfar1222dotx(complex<double> x[], complex<double> y[]);
void mvm_BDPmatrix(std::complex<double> x[], std::complex<double> y[]);

void stage_solution(){

    int itermax=val_itermax;
    int iterrestart=val_itermax;
    double soltolerance=val_soltol;
    
    double delta_theta=(source_PW[1]-source_PW[0]);
    double delta_phi=(source_PW[3]-source_PW[2]);
    complex <double> val_RHS_ratio;
    complex <double> *vec_RHS_prev=new complex <double>[num_rwgs*2];
    dir_PW = new double[4];
    aggregation=new complex<double>[size_agg*2];
    disaggregation=new complex<double>[size_agg*2];
    vec_sol=new complex<double>[num_rwgs*2];
    for(int ind0=0;ind0<(num_rwgs*2);ind0++){
        vec_sol[ind0]=complex<double>(0.0, 0.0);
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    int ind_sol=0;
    for(int ind_theta=0; ind_theta<num_source[0]; ind_theta++){
        for(int ind_phi=0; ind_phi<num_source[1]; ind_phi++){
        calculate_dir_pw(num_source, delta_theta, delta_phi, dir_PW, ind_theta, ind_phi);
        calculate_rhsvector(dir_PW);
        if (val_pc==1){mvm_BDPmatrix(vec_RHS, vec_RHS);}
        calculate_initguess(ind_sol, val_RHS_ratio, vec_RHS_prev);

        if((val_pc==0)||(val_pc==1)) gmres((num_rwgs*2),vec_RHS,itermax,soltolerance,vec_sol,0);
        if(val_pc==2) fgmres((num_rwgs*2),vec_RHS,itermax,iterrestart,soltolerance,vec_sol);

        copy_rhsvector(vec_RHS_prev);
        calculate_farfield(ind_sol);
        ind_sol++;
        }
    }



    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    printf("Solution Part Takes %f s \n", time_span.count());

    std::cout.flush();

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
        for(int ind_rwg=0;ind_rwg<(num_rwgs*2);ind_rwg++){
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
    for(int ind_rwg=0;ind_rwg<(num_rwgs*2);ind_rwg++){
        vec_RHS_prev[ind_rwg]=vec_RHS[ind_rwg];
    }
}