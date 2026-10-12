#include "matrixBlock.h"
#include "matrixoperations.h"
#include "matrix_inverse.h"
#include <iostream>

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
extern int num_filledboxesall;
extern int *num_children;
extern matrixBlock* blockdiagonalmatrix11;
extern matrixBlock* blockdiagonalmatrix12;
extern matrixBlock* blockdiagonalmatrix21;
extern matrixBlock* blockdiagonalmatrix22;

void calculate_blockdiagonalpreconditioner(){
    complex<double>** mattest;
    complex<double>** result;
    int matsize;

    printf("num_filledboxesall: %d \n", num_filledboxesall);

    blockdiagonalmatrix11 = new matrixBlock[num_filledboxes[num_mlfmalevels-1]];
    blockdiagonalmatrix12 = new matrixBlock[num_filledboxes[num_mlfmalevels-1]];
    blockdiagonalmatrix21 = new matrixBlock[num_filledboxes[num_mlfmalevels-1]];
    blockdiagonalmatrix22 = new matrixBlock[num_filledboxes[num_mlfmalevels-1]];

    int box_start=num_filledboxes[num_mlfmalevels*2-1];
    int box_end=num_filledboxes[num_mlfmalevels*3-1]+1;

    for(int box_test=box_start;box_test<box_end;box_test++){
        for(int ind_nearbox=0;ind_nearbox<27;ind_nearbox++){
            if(boxes_near[box_test*27+ind_nearbox]!=-1){
                int box_basis=boxes_near[box_test*27+ind_nearbox];
                if(box_test==box_basis){
                    int num_rwgs_inbox = num_children[box_test];
                    blockdiagonalmatrix11[box_test-box_start].initmatrixBlock(num_rwgs_inbox,num_rwgs_inbox);
                    blockdiagonalmatrix12[box_test-box_start].initmatrixBlock(num_rwgs_inbox,num_rwgs_inbox);
                    blockdiagonalmatrix21[box_test-box_start].initmatrixBlock(num_rwgs_inbox,num_rwgs_inbox);
                    blockdiagonalmatrix22[box_test-box_start].initmatrixBlock(num_rwgs_inbox,num_rwgs_inbox);
                    int rwg_test_start=children[box_test*2+0];
                    int rwg_test_end=children[box_test*2+1]+1;
                    long int ind_Znear=nearindex[(box_test-box_start)*27+ind_nearbox];
                    for(int ind_rwg_test=rwg_test_start;ind_rwg_test<rwg_test_end;ind_rwg_test++){
                        int testfunction=rwg_ind[ind_rwg_test];
                        blockdiagonalmatrix11[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction; 
                        blockdiagonalmatrix11[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction; //Test and basis functions are the same for self-box
                        blockdiagonalmatrix12[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction; 
                        blockdiagonalmatrix12[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction; //Test and basis functions are the same for self-box
                        blockdiagonalmatrix21[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction; 
                        blockdiagonalmatrix21[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction; //Test and basis functions are the same for self-box
                        blockdiagonalmatrix22[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction; 
                        blockdiagonalmatrix22[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction; //Test and basis functions are the same for self-box
                        /*
                        blockdiagonalmatrix11[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction; 
                        blockdiagonalmatrix11[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction; //Test and basis functions are the same for self-box
                        blockdiagonalmatrix12[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction; 
                        blockdiagonalmatrix12[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction+num_rwgs; //Test and basis functions are the same for self-box
                        blockdiagonalmatrix21[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction+num_rwgs; 
                        blockdiagonalmatrix21[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction; //Test and basis functions are the same for self-box
                        blockdiagonalmatrix22[box_test-box_start].indrows[ind_rwg_test-rwg_test_start]=testfunction+num_rwgs; 
                        blockdiagonalmatrix22[box_test-box_start].indcols[ind_rwg_test-rwg_test_start]=testfunction+num_rwgs; //Test and basis functions are the same for self-box
                        */
                    }
                    for(int ind_rwg_test=rwg_test_start;ind_rwg_test<rwg_test_end;ind_rwg_test++){
                        int testfunction=rwg_ind[ind_rwg_test];
                        int rwg_basis_start=children[box_basis*2+0];
                        int rwg_basis_end=children[box_basis*2+1]+1;
                        for(int ind_rwg_basis=rwg_basis_start;ind_rwg_basis<rwg_basis_end;ind_rwg_basis++){
                            int basisfunction=rwg_ind[ind_rwg_basis];
                            blockdiagonalmatrix11[box_test-box_start].assignvalue(ind_rwg_test-rwg_test_start,ind_rwg_basis-rwg_basis_start,matrix_Znear11[ind_Znear]);
                            blockdiagonalmatrix12[box_test-box_start].assignvalue(ind_rwg_test-rwg_test_start,ind_rwg_basis-rwg_basis_start,matrix_Znear12[ind_Znear]);
                            blockdiagonalmatrix21[box_test-box_start].assignvalue(ind_rwg_test-rwg_test_start,ind_rwg_basis-rwg_basis_start,matrix_Znear21[ind_Znear]);
                            blockdiagonalmatrix22[box_test-box_start].assignvalue(ind_rwg_test-rwg_test_start,ind_rwg_basis-rwg_basis_start,matrix_Znear22[ind_Znear]);
                            ind_Znear++;
                        }
                    }
                }
            }
        }
    }



    printf("copying finished \n");

    for(int box_test=box_start;box_test<box_end;box_test++){
        mattest = blockdiagonalmatrix11[box_test-box_start].submatrix;  
        matsize = blockdiagonalmatrix11[box_test-box_start].sizerow;
        complex<double>** result = matrix_inverse(mattest, matsize);
        for(int indr=0; indr<matsize; indr++){
            for(int indc=0; indc<matsize; indc++){
                blockdiagonalmatrix11[box_test-box_start].assignvalue(indr,indc,result[indr][indc]);
            }
        }
        for(int indr=0; indr<matsize; indr++)
            delete[] result[indr];
    }
    for(int box_test=box_start;box_test<box_end;box_test++){
        mattest = blockdiagonalmatrix12[box_test-box_start].submatrix;  
        matsize = blockdiagonalmatrix12[box_test-box_start].sizerow;
        complex<double>** result = matrix_inverse(mattest, matsize);
        for(int indr=0; indr<matsize; indr++){
            for(int indc=0; indc<matsize; indc++){
                blockdiagonalmatrix12[box_test-box_start].assignvalue(indr,indc,result[indr][indc]);
            }
        }
        for(int indr=0; indr<matsize; indr++)
            delete[] result[indr];
    }
    for(int box_test=box_start;box_test<box_end;box_test++){
        mattest = blockdiagonalmatrix21[box_test-box_start].submatrix;  
        matsize = blockdiagonalmatrix21[box_test-box_start].sizerow;
        complex<double>** result = matrix_inverse(mattest, matsize);
        for(int indr=0; indr<matsize; indr++){
            for(int indc=0; indc<matsize; indc++){
                blockdiagonalmatrix21[box_test-box_start].assignvalue(indr,indc,result[indr][indc]);
            }
        }
        for(int indr=0; indr<matsize; indr++)
            delete[] result[indr];
    }
    for(int box_test=box_start;box_test<box_end;box_test++){
        mattest = blockdiagonalmatrix22[box_test-box_start].submatrix;  
        matsize = blockdiagonalmatrix22[box_test-box_start].sizerow;
        complex<double>** result = matrix_inverse(mattest, matsize);
        for(int indr=0; indr<matsize; indr++){
            for(int indc=0; indc<matsize; indc++){
                blockdiagonalmatrix22[box_test-box_start].assignvalue(indr,indc,result[indr][indc]);
            }
        }
        for(int indr=0; indr<matsize; indr++)
            delete[] result[indr];
    }





/*
    int num_boxes=num_filledboxes[num_mlfmalevels-1];
    for(int box_test=0;box_test<10;box_test++){
        int num_test = blockdiagonalmatrix21[box_test].sizerow;
        int num_basis = blockdiagonalmatrix21[box_test].sizecol;
        printf("%d %d\n", num_test, num_basis);
        for(int i=0; i<num_test; i++){
            for(int j=0; j<num_basis; j++){
                complex<double> matval = blockdiagonalmatrix21[box_test].submatrix[i][j];
                printf("%d %d %1.8f %1.8f \n", i, j, real(matval), imag(matval));
            }
        }
    }
*/



    printf("inverse operation finished \n");

    //delete[] blockdiagonalmatrix;

    std::cout.flush();
    
}