#include <complex>

using namespace std;

extern int num_rwgs;
extern int val_pc;

void mvm_Zneardotx(complex<double> x[], complex<double> y[]);
void mvm_Zfar1121dotx(complex<double> x[], complex<double> y[]);
void mvm_Zfar1222dotx(complex<double> x[], complex<double> y[]);
void mvm_BDPmatrix(std::complex<double> x[], std::complex<double> y[]);

void mvm(complex<double> x[], complex<double> y[]){

    mvm_Zneardotx(x,y);
    mvm_Zfar1121dotx(x,y);
    mvm_Zfar1222dotx(x,y);

    if ((val_pc==1)){mvm_BDPmatrix(y, y);}

/*
    if ((val_pc==4)){
        complex <double> *y_temp=new complex <double>[num_rwgs*2];
        for(int i=0; i<(num_rwgs*2); i++)
            y_temp[i]=y[i];
        for(int i=0; i<(num_rwgs*2); i++)
            y[i]=complex<double>(0.0, 0.0);
        mvm_Zneardotx(y_temp,y);
        mvm_Zfar1121dotx(y_temp,y);
        mvm_Zfar1222dotx(y_temp,y);
    }
*/
    
}
