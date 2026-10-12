#include <complex>

using namespace std;

extern int num_rwgs;
extern int val_pc;

void mvm_Zneardotx(complex<double> x[], complex<double> y[]);
void mvm_Zfar1121dotx(complex<double> x[], complex<double> y[]);
void mvm_Zfar1222dotx(complex<double> x[], complex<double> y[]);
void mvm_BDPmatrix(std::complex<double> x[], std::complex<double> y[]);

void mvm_physical(complex<double> x[], complex<double> y[]){

    mvm_Zneardotx(x,y);
    mvm_Zfar1121dotx(x,y);
    mvm_Zfar1222dotx(x,y);

    
}

// Left preconditioning is excluded when checking the physical equation.
void mvm(complex<double> x[], complex<double> y[]){
    mvm_physical(x, y);
    if (val_pc==1){mvm_BDPmatrix(y, y);}
}
