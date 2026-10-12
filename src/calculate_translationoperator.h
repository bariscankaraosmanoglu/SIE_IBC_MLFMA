#include <cmath>
#include <complex>
#include <iostream>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

extern complex<double> *trans_op; // created in here
extern complex<double> *trans_op_pc; // created in here
extern double *size_boxes;
extern double val_kO;
extern int *num_truncation;
extern int *num_truncation_pc;
extern int num_mlfmalevels;
extern int *trans_index;
extern int *trans_index_pc;
extern int size_trans;
extern int size_trans_pc;
extern double val_approx;
extern int val_pc;
extern int parallel_omp_on;
extern int num_threads;
extern FILE *file_timer;

void phi_locations_weights(int , double*, double*);
void theta_locations_weights(int , double*, double*);
void calculate_sin_cos(int, int, double*, double*, double*, double*, double*, double*);
double hankel_func(double, int);
