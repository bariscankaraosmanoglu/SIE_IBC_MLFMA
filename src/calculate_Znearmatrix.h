#include "solver_globals.h"
#include <stdio.h>
#include <cmath>
#include <complex>
#include <chrono>
#include <omp.h>
#include <iostream>

using namespace std;
using namespace std::chrono;

extern double val_etaO;
extern double val_cfiealpha[];
extern complex<double> val_impedance[];
extern int *children;
extern double val_kO;
extern double val_epsilonO;
extern double *vec_trinormals;
extern int *smallestbox2tri;
extern int val_maxtri;
extern double pi;
extern double val_muO;
extern int *boxes_near;
extern complex<double> *matrix_Znear11; // created in here
extern complex<double> *matrix_Znear12; // created in here
extern complex<double> *matrix_Znear21; // created in here
extern complex<double> *matrix_Znear22; // created in here
extern long int *nearindex;
extern double *node_pos;
extern int *num_children;
extern int num_mlfmalevels;
extern long int size_Znear;
extern int *num_filledboxes;
extern double val_omega;
extern int *rwg2box;
extern double *rwg_lengths;
extern int *rwg_reverseind;
extern double *tri_areas;
extern int *tri_nodes;
extern int *tri_media;
extern int *tri2rwgindex;
extern int *tri2rwgnode;
extern int *tri2rwgsign;
extern int parallel_omp_on;
extern int num_threads;
extern FILE *file_timer;


void quadrature_points(double r1[],double r2[],double r3[],double weights[],double locations[][3],int num_points);
void tri_analytical_int(double P1[],double P2[],double P3[],double Pobs[],double triarea,double intval1[],double intval2[],double intval3[],double intval4[]);
double dotProduct3d(double vec1[], double vec2[]);
void subtract3d(double vec1[], double vec2[], double vec3[]);
