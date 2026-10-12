#include <stdio.h>
#include <stdlib.h>
#include <cmath>

extern int *boxes_near;
extern long int *nearindex; // created here
extern int *num_children;
extern int num_mlfmalevels;
extern int num_rwgs;
extern int *num_filledboxes;
extern int *rwg2box; // created here
extern int *rwg_ind;
extern int *rwg_reverseind; // created here
extern int *rwg_tris;
extern long int size_Znear; // created here
extern int *smallestbox2tri; // created here
extern int val_maxtri; // created here

int compareint (const void *a, const void *b);
