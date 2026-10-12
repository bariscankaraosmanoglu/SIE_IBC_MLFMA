#include "solver_globals.h"
#include "calculate_loadbalance.h"
#include <complex>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

extern int *children;
extern int *boxes_near;
extern complex<double> *matrix_Znear11;
extern complex<double> *matrix_Znear12;
extern complex<double> *matrix_Znear21;
extern complex<double> *matrix_Znear22;
extern long int *nearindex;
extern int num_mlfmalevels;
extern int *num_filledboxes;
extern int num_rwgs;
extern int *rwg_ind;
extern int parallel_omp_on;
extern int num_threads;
extern FILE *file_timer;

extern long int size_Znear;

void mvm_Zneardotx(complex<double> x[], complex<double> y[]){
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    int box_start=num_filledboxes[num_mlfmalevels*2-1];
    int box_end=num_filledboxes[num_mlfmalevels*3-1]+1;

//    long int ind_Znear=0;

    omp_set_num_threads(num_threads);
    #pragma omp parallel if(parallel_omp_on)
    {

        const auto box_range=near_box_range(omp_get_thread_num(),omp_get_num_threads());
        double timer_start = omp_get_wtime();
        double timer_total = 0.0;


//        #pragma omp for schedule(dynamic,50)
//        for(int box_test=box_start;box_test<box_end;box_test++){
        for(int box_test=box_range.first;box_test<box_range.second;box_test++){
            double timer_in = omp_get_wtime();
            for(int ind_nearbox=0;ind_nearbox<27;ind_nearbox++){
                if(boxes_near[box_test*27+ind_nearbox]!=-1){
                    int box_basis=boxes_near[box_test*27+ind_nearbox];
                    int rwg_test_start=children[box_test*2+0];
                    int rwg_test_end=children[box_test*2+1]+1;
                    long int ind_Znear=nearindex[(box_test-box_start)*27+ind_nearbox];
                    for(int ind_rwg_test=rwg_test_start;ind_rwg_test<rwg_test_end;ind_rwg_test++){
                        int testfunction=rwg_ind[ind_rwg_test];
                        int rwg_basis_start=children[box_basis*2+0];
                        int rwg_basis_end=children[box_basis*2+1]+1;
                        for(int ind_rwg_basis=rwg_basis_start;ind_rwg_basis<rwg_basis_end;ind_rwg_basis++){
                            int basisfunction=rwg_ind[ind_rwg_basis];
                            y[testfunction]+=matrix_Znear11[ind_Znear]*x[basisfunction];
                            y[testfunction]+=matrix_Znear12[ind_Znear]*x[basisfunction+num_rwgs];
                            y[testfunction+num_rwgs]+=matrix_Znear21[ind_Znear]*x[basisfunction];
                            y[testfunction+num_rwgs]+=matrix_Znear22[ind_Znear]*x[basisfunction+num_rwgs];
                            ind_Znear++;
                        }
                    }
                }
            }
            timer_total+=omp_get_wtime()-timer_in;
        }

        double timer_end = omp_get_wtime();
//        fprintf(file_timer,"MVMNear %d %1.5f %1.5f\n", omp_get_thread_num(), timer_end-timer_start, timer_total );

    }



    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
//    printf("Znear MVM Calculation Takes %f s \n", time_span.count());


}



