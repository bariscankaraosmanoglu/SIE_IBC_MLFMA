#ifndef SOLVER_GLOBALS_H
#define SOLVER_GLOBALS_H
#include <complex>
// Shared storage types must agree across assembly, operators, and cleanup.
extern char modelname[200], modeltriname[200], modelnodename[200];
extern double val_cfiealpha[50], source_PW[4];
extern std::complex<double> val_impedance[50];
extern long int *nearindex, *agg_index, *agg_index_pc;
#endif
