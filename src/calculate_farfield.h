#include <stdio.h>
#include <cmath>
#include <complex>
#include <chrono>
#include <iostream>

using namespace std;
using namespace std::chrono;

extern int num_nodes;
extern int num_tris;
extern double *node_pos;
extern int *tri_nodes;

extern int *num_trispernode;
extern int *nodetrisindex;
extern int num_rwgs;
extern int *rwg_tris;
extern int *rwg_nodes;
extern double *rwg_lengths;
extern int *tri2rwgindex;
extern int *tri2rwgnode;
extern int *tri2rwgsign;
extern double *tri_areas;

extern double pi;
extern double val_e;
extern double val_epsilonO;
extern double val_muO;
extern double val_freq;
extern double val_omega;
extern double val_c;
extern double val_lambdaO;
extern double val_kO;
extern double val_etaO;
extern double *dir_PW;
extern int val_typeff;
extern int *num_ffdir_tp;
extern double dir_ff;
extern int num_ffdir_cut;
extern int num_ffdir;
extern double *ffdir_phi;
extern double *ffdir_theta;

extern complex<double> *vec_sol;

extern complex<double> *vec_farfield;

void quadrature_points(double r1[],double r2[],double r3[],double weights[],double locations[][3],int num_points);