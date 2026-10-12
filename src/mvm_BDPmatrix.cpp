#include <complex>
#include "matrixBlock.h"

extern matrixBlock* blockdiagonalmatrix11;
extern matrixBlock* blockdiagonalmatrix12;
extern matrixBlock* blockdiagonalmatrix21;
extern matrixBlock* blockdiagonalmatrix22;
extern int num_rwgs;
extern int *num_filledboxes;
extern int num_mlfmalevels;

extern int *children;
extern int *boxes_near;
extern long int *nearindex;
extern int num_mlfmalevels;
extern int *num_filledboxes;
extern int num_rwgs;
extern int *rwg_ind;

void mvm_BDPmatrix(std::complex<double> x[], std::complex<double> y[]){

    std::complex<double>* y2 = new std::complex<double>[num_rwgs*2];
    for(int i=0; i<(num_rwgs*2); i++)
        y2[i]=complex<double> (0.0, 0.0);
    
    // Apply Block-Diagonal Preconditioner


    int num_boxes=num_filledboxes[num_mlfmalevels-1];


    for(int box_test=0;box_test<num_boxes;box_test++){
        int num_test = blockdiagonalmatrix11[box_test].sizerow;
        int num_basis = blockdiagonalmatrix11[box_test].sizecol;
        for(int i=0; i<num_test; i++){
            int indr=blockdiagonalmatrix11[box_test].indrows[i];
            for(int j=0; j<num_basis; j++){
                int indc=blockdiagonalmatrix11[box_test].indcols[j];
                y2[indr]+=blockdiagonalmatrix11[box_test].submatrix[i][j]*x[indc];
            }
        }
    }


/*
    for(int box_test=0;box_test<num_boxes;box_test++){
        int num_test = blockdiagonalmatrix12[box_test].sizerow;
        int num_basis = blockdiagonalmatrix12[box_test].sizecol;
        for(int i=0; i<num_test; i++){
            int indr=blockdiagonalmatrix12[box_test].indrows[i];
            for(int j=0; j<num_basis; j++){
                int indc=blockdiagonalmatrix12[box_test].indcols[j]+num_rwgs;
                y2[indr]+=blockdiagonalmatrix12[box_test].submatrix[i][j]*x[indc];
            }
        }
    }
*/

/*
    for(int box_test=0;box_test<num_boxes;box_test++){
        int num_test = blockdiagonalmatrix21[box_test].sizerow;
        int num_basis = blockdiagonalmatrix21[box_test].sizecol;
        for(int i=0; i<num_test; i++){
            int indr=blockdiagonalmatrix21[box_test].indrows[i]+num_rwgs;
            for(int j=0; j<num_basis; j++){
                int indc=blockdiagonalmatrix21[box_test].indcols[j];
                y2[indr]+=blockdiagonalmatrix21[box_test].submatrix[i][j]*x[indc];
            }
        }
    }
*/

    for(int box_test=0;box_test<num_boxes;box_test++){
        int num_test = blockdiagonalmatrix22[box_test].sizerow;
        int num_basis = blockdiagonalmatrix22[box_test].sizecol;
        for(int i=0; i<num_test; i++){
            int indr=blockdiagonalmatrix22[box_test].indrows[i]+num_rwgs;
            for(int j=0; j<num_basis; j++){
                int indc=blockdiagonalmatrix22[box_test].indcols[j]+num_rwgs;
                y2[indr]+=blockdiagonalmatrix22[box_test].submatrix[i][j]*x[indc];
            }
        }
    }


    for(int i=0; i<(num_rwgs*2); i++)
        y[i]=y2[i];

    delete[] y2;
    

}