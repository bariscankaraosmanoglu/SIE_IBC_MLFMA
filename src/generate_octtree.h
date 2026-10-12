#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

using namespace std;

extern int num_nodes;
extern int num_tris;
extern int num_rwgs;
extern double *node_pos;
extern int *tri_nodes;
extern int *rwg_nodes;
extern double val_lambdaO;
extern double val_kO;

extern int val_pc;
extern double val_approx;

extern int num_mlfmalevels;
extern double *size_boxes;
extern int *num_truncation;
extern int *num_truncation_pc;

extern int num_digits;

extern int num_filledboxesall;
extern int *num_filledboxes;
extern int *num_children;
extern double *box_pos;
extern int *children;
extern int *children2parent;
extern int *rwg_ind;

struct sortedElements {
    int elemID;
    long int elemVal;
};

int compare (const void *a, const void *b);
int compareind (const void *a, const void *b);
long int dec2bin2dec(int,int,int,int);
