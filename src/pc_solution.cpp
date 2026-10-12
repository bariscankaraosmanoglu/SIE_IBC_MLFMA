#include <stdio.h>
#include <complex>
#include <chrono>

using namespace std;
using namespace std::chrono;

extern int num_rwgs;
extern int val_pc;
extern int val_itermax_pc;

int fgmres_in(int size_vecRHS,complex<double> vec_RHS[],int itermax,int iterrestart,double soltolerance,complex<double> vec_sol[]);
int gmres(int size_vecRHS,complex<double> vec_RHS_i[],int itermax,double soltolerance,complex<double> vec_sol_i[], int pc_in);

void pc_solution(complex<double> vec_RHS_i[], complex<double> vec_sol_i[], int size_vecRHS){


    int itermax=val_itermax_pc;
    int iterrestart=val_itermax_pc;
    double soltolerance=0.5;

    for(int ind0=0;ind0<size_vecRHS;ind0++){
        vec_sol_i[ind0]=complex<double>(0.0, 0.0);
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    fgmres_in(size_vecRHS,vec_RHS_i,itermax,iterrestart,soltolerance,vec_sol_i);
    //stat_solution++; // for nothing

    //gmres(size_vecRHS,vec_RHS_i,itermax,soltolerance,vec_sol_i,1);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
//    printf("Preconditioner Part Takes %f s \n", time_span.count());



}



//#include <stdio.h>
//#include <complex>
//
//using namespace std;
//
//void pc_solution(complex<double> vec_RHS_i[], complex<double> vec_sol_i[], int size_vecRHS){
//
//    for(int ind=0;ind<size_vecRHS;ind++){
//        vec_sol_i[ind]=vec_RHS_i[ind];
//    }
//
//}

