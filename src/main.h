#include "solver_globals.h"
#include "calculate_loadbalance.h"
#include <stdio.h>
#include <complex>
#include <chrono>
#include "matrixBlock.h"

using namespace std;
using namespace std::chrono;

// Variables
// created in read_simulationinput
//const char *modeltriname;
//const char *modelnodename;
char modeltriname[200];
char modelnodename[200];
char modelname[200];
int num_mlfmalevels;
int num_digits;
int num_interpol;
int source_type;
int val_itermax;
int val_pc;
int val_itermax_pc;
int num_media;
complex<double> val_impedance[50];
double val_cfiealpha[50];
double val_soltol;
double val_approx;
double dir_cfc[2];
complex<double> *source_CSP;
complex<double> *source_current;
complex<double> *source_DG;
double *source_dipole;
double source_PW[4];
double *dir_PW;
// Source counts have static ownership; num_source is a borrowed view.
int num_source_storage[2]{};
int *num_source=num_source_storage;
int val_pol;
int parallel_omp_on;
int num_threads;
FILE *file_timer;
double dir_ff;
int num_freqs;
double freq_start, freq_end;

// created in read_meshmodel
int num_nodes=0;
int num_tris=0;
double *node_pos;
int *tri_nodes;
int *tri_media;

// created in generate_edgeconnection
int *nodetrisindex;
int num_rwgs=0;
int *rwg_tris;
int *rwg_nodes;
double *rwg_lengths;
int *tri2rwgindex;
int *tri2rwgnode;
int *tri2rwgsign;
double *tri_areas;
double *vec_trinormals;

// created in generate_octtree
int *num_truncation;
int *num_truncation_pc;
double *size_boxes;
int num_filledboxesall;
int *num_filledboxes;
int *num_children;
double *box_pos;
int *children;
int *children2parent;
int *rwg_ind;

// created in genetare_constants
double pi=3.1415926535897932384626433;
double val_e=2.7182818284590452353602874;
double val_epsilonO=8.854187817620390e-12;
double val_muO=1.2566370621219e-6;
double val_freq;
double val_omega;
double val_c;
double val_lambdaO;
double val_kO;
double val_etaO;
int val_typeff;
int *num_ffdir_tp;
int num_ffdir_cut;
int num_ffdir;
double *ffdir_phi;
double *ffdir_theta;

// created in generate_nearfarboxes
int *boxes_far;
int *boxes_near;

// created in calculate_loadbalance

// created in generate_Znearcomponent
long int *nearindex;
int *rwg2box;
int *rwg_reverseind;
long int size_Znear;
int *smallestbox2tri;
int val_maxtri;

// created in generate_Zfarcomponents
long int *agg_index;
long int *agg_index_pc;
long int size_agg;
long int size_agg_pc;

// created in generate_translationcomponent
int *trans_index;
int *trans_index_pc;
int size_trans;
int size_trans_pc;

// created in calculate_translationoperator
complex<double> *trans_op;
complex<double> *trans_op_pc;

// created in calculate_Znearmatrix
complex<double> *matrix_Znear11;
complex<double> *matrix_Znear12;
complex<double> *matrix_Znear21;
complex<double> *matrix_Znear22;

// created in calculate_blockdiagonalpreconditioner
matrixBlock* blockdiagonalmatrix11 = nullptr;
matrixBlock* blockdiagonalmatrix22 = nullptr;

// created in calculate_rhsvector
complex<double> *vec_RHS;

// created in
complex<double> *aggregation;
complex<double> *disaggregation;

// created in stage_solution
complex<double> *vec_sol;

// created in calculate_farfield
complex<double> *vec_farfield;


// Functions
void read_simulationinput();
void read_meshmodel();
void generate_edgeconnection();
void genetare_constants();
void generate_octtree();
void generate_nearfarboxes();
void generate_Znearcomponent();
void generate_Zfarcomponents();
void generate_translationcomponent();
void calculate_loadbalance();
void calculate_translationoperator();
void calculate_Znearmatrix();
void calculate_blockdiagonalpreconditioner();
// void calculate_rhsvector();

void mvm_Zneardotx();
void mvm_Zfardotx();

int stage_solution();
void stage_aftersolution();

void clear_memory();
