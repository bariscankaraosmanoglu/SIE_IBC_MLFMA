#include <complex>

using namespace std;

extern int num_rwgs;
extern int val_pc;
extern double val_approx;

void mvm_Zneardotx(complex<double> x[], complex<double> y[]);
void mvm_Zfar1121dotx_in(complex<double> x[], complex<double> y[]);
void mvm_Zfar1222dotx_in(complex<double> x[], complex<double> y[]);

void mvm_in(complex<double> x[], complex<double> y[]){

    mvm_Zneardotx(x,y);
    if(val_approx>0){
        mvm_Zfar1121dotx_in(x,y);
        mvm_Zfar1222dotx_in(x,y);
    }

}
