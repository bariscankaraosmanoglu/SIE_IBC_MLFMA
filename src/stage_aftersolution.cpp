#include <stdio.h>
#include <complex>

using namespace std;

void writeout_rwgcoefficients();
void calculate_farfield();
void calculate_currentontris();
void calculate_currentfieldcontribution();

extern int num_rwgs;
extern complex<double> *vec_sol;

void stage_aftersolution(){

    writeout_rwgcoefficients();
    // calculate_farfield();
    calculate_currentontris();
    calculate_currentfieldcontribution();

}

void writeout_rwgcoefficients(){
    FILE *fptr;
    fptr=fopen("result_rwgcoefficient.txt","w");
    for(int ind=0;ind<(num_rwgs*2);ind++){
        fprintf(fptr,"%1.15f %1.15f \n", vec_sol[ind].real(), vec_sol[ind].imag());
    }
    fclose(fptr);
}
