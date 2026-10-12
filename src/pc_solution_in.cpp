#include <stdio.h>
#include <complex>
#include <chrono>

using namespace std;
using namespace std::chrono;

extern int num_rwgs;

int fgmres_in_in(int size_vecRHS,complex<double> vec_RHS[],int itermax,int iterrestart,double soltolerance,complex<double> vec_sol[]);
int gmres(int size_vecRHS,complex<double> vec_RHS_i[],int itermax,double soltolerance,complex<double> vec_sol_i[], int pc_in);

void pc_solution_in(complex<double> vec_RHS_i[], complex<double> vec_sol_i[], int size_vecRHS){

    
    int itermax=20;
    int iterrestart=20;
    double soltolerance=0.4;

    for(int ind0=0;ind0<num_rwgs;ind0++){
        vec_sol_i[ind0]=complex<double>(0.0, 0.0);
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    gmres(size_vecRHS,vec_RHS_i,itermax,soltolerance,vec_sol_i,2);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    

/*
    for(int ind=0;ind<size_vecRHS;ind++){
        vec_sol_i[ind]=vec_RHS_i[ind];
    }
*/

}

