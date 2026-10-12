#include <stdio.h>
#include <complex>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <unistd.h>

using namespace std;

extern int parallel_omp_on;

void mvm(complex<double> x[], complex<double> y[]);
void mvm_in(complex<double> x[], complex<double> y[]);
void mvm_Zneardotx(complex<double> x[], complex<double> y[]);
double LANORM(complex<double> vec_in[],int size_vec);
void cleartozero(complex<double> vec_in[],int size_vec);
void destructor2D(complex<double> **data_in, int num_row, int num_col);
complex<double> ** constructor2D(int num_row, int num_col);
complex<double> SCALARSIGN(complex<double> inputval);
complex<double> SCALARP(complex<double> vec_in1[],complex<double> vec_in2[],int size_vec);
complex<double>* SOLVELSP(complex<double>** A, complex<double>* b, int sizeA);
//void SOLVELSP(complex<double>** A, complex<double>* b, int sizeA, complex<double>* x);

int gmres(int size_vecRHS,complex<double> vec_RHS[],int itermax,double soltolerance,complex<double> vec_sol[], int pc_in){

    int n = size_vecRHS;
    int flag;
    int imin=0;
    int jmin=0;
    int evalxm=0;
    int stag=0;
    int moresteps=0;
    int maxmsteps=itermax;
    int maxstagsteps=3;
    int minupdated=0;
    double relres;
    double resvec[itermax+1];
    double normr=0.0;
    double normr_act=0.0;
    double normrmin;
    double normU=0.0;
    double normv=0.0;
    double normUhouse=0.0;
    double rho=0.0;
    complex<double> alpha;
    complex<double>** Jgrot;
    complex<double>** Ugrot;
    complex<double>** Rgrot;
    complex<double>* wgrot;
    complex<double>* Uhouse=new complex<double>[n];
    complex<double>* v=new complex<double>[n];
    complex<double>* vec_mvm=new complex<double>[n];
    complex<double>* r=new complex<double>[n];
    complex<double>* additive=new complex<double>[n];

    complex<double> beta;
    complex<double> dummyconj;
    complex<double> tmpv;
    complex<double> val_product;
    complex<double>*vec_temp;
    //complex<double>* ytmp=new complex<double> [itermax];
    //cleartozero(vec_mvm,n);
    //cleartozero(r,n);


    int maxit = itermax;
    if (itermax>n) maxit = n;
    int inner= maxit;

    // Check for all zero right hand side vector => all zero solution
    double n2b = LANORM(vec_RHS,n);
    if(n2b==0.0){
        for(int ind1=0; ind1<n; ind1++){
            vec_sol[ind1]=complex<double>(0.0, 0.0);
        }
        flag=0;
        relres=0.0;
        resvec[0]=0.0;
        return 1;
    }
    double tolb = soltolerance*n2b;

    //mvm(vec_sol,vec_mvm);
    for(int i=0; i<n; i++)
        vec_mvm[i]=complex<double>(0.0, 0.0);
        
    if (pc_in==0){mvm(vec_sol,vec_mvm);}
    if (pc_in==1){mvm_in(vec_sol,vec_mvm);}
    if (pc_in==2){mvm_Zneardotx(vec_sol,vec_mvm);}

    for(int i=0; i<n; i++)
        r[i]=vec_RHS[i]-vec_mvm[i];

    normr = LANORM(r,n);

    if(normr<=tolb){
        //printf("GMRES LOCATION FOR THIS ERROR: 1\n");
        return 0;
    }
    resvec[0]=normr;
    normrmin=normr;

    //Allocate J to hold the Given's rotation constants
    Jgrot=constructor2D(2, inner);
    Ugrot=constructor2D(inner+1, n); // switched the dimensions
    Rgrot=constructor2D(inner, inner);
    wgrot=new complex<double>[inner+1];
    //cleartozero(wgrot,(inner+1));
    
    for(int i=0; i<n; i++)
        Uhouse[i]=r[i];
    
    beta = SCALARSIGN(r[0])*normr; 
    Uhouse[0]+=beta;
    normU = LANORM(Uhouse,n);

    for(int i=0; i<n; i++){
        Uhouse[i]/=normU;
        Ugrot[0][i]=Uhouse[i];
    }
    wgrot[0] = -beta;

    for(int initer=0; initer<inner; initer++){
        dummyconj=conj(Uhouse[initer]);
        for(int ind1=0; ind1<n; ind1++)
            v[ind1]=-2.0*dummyconj*Uhouse[ind1];
        v[initer]+=complex<double>(1.0, 0.0);

        for(int k=(initer-1); k>=0; k--){
            vec_temp=&Ugrot[k][0];
            val_product=2.0*SCALARP(vec_temp,v,n);
            for(int ind1=0; ind1<n; ind1++){
                v[ind1]-=val_product*Ugrot[k][ind1];
            }
        }
        normv = LANORM(v,n);

        for(int ind1=0; ind1<n; ind1++){
            v[ind1]/=normv;
        }

        //cleartozero(vec_mvm,n);
        //mvm(v,vec_mvm);
        for(int i=0; i<n; i++)
            vec_mvm[i]=complex<double>(0.0, 0.0);
            
        if (pc_in==0){mvm(v,vec_mvm);}
        if (pc_in==1){mvm_in(v,vec_mvm);}
        if (pc_in==2){mvm_Zneardotx(v,vec_mvm);}

        for(int ind1=0; ind1<n; ind1++){
            v[ind1]=vec_mvm[ind1];
        }

        // Form Pj*Pj-1*...P1*Av.
        for(int k=0; k<=initer; k++){
            vec_temp=&Ugrot[k][0];
            val_product=2.0*SCALARP(vec_temp,v,n);
            for(int ind1=0; ind1<n; ind1++){
                v[ind1]-=val_product*Ugrot[k][ind1];
            }
        }

        //Determine Pj+1
        if(initer!=n){
            for(int ind1=0; ind1<n; ind1++)
                Uhouse[ind1]=v[ind1];
            for(int ind1=0; ind1<=initer; ind1++)
                Uhouse[ind1]=complex<double>(0.0, 0.0);
            normUhouse=LANORM(Uhouse,n);
            
            if(normUhouse!=0.0){
                alpha=SCALARSIGN(v[initer+1])*normUhouse;
                Uhouse[initer+1]=Uhouse[initer+1]+alpha;
                normU=LANORM(Uhouse,n);

                for(int i=0;i<n;i++)
                    Uhouse[i]/=normU;
                for(int i=0;i<n;i++)
                    Ugrot[initer+1][i]=Uhouse[i];
                for(int i=initer+2;i<n;i++)
                    v[i]=complex<double>(0.0, 0.0);
                v[initer+1]=-alpha;
            }
        }

        //Apply Given's rotations to the newly formed v.
        for(int colJ=0; colJ<=(initer-1); colJ++){
            tmpv=v[colJ];
            v[colJ]=conj(Jgrot[0][colJ])*tmpv+conj(Jgrot[1][colJ])*v[colJ+1];
            v[colJ+1]=-Jgrot[1][colJ]*tmpv+Jgrot[0][colJ]*v[colJ+1];
        }

        //Compute Given's rotation Jm
        if(initer!=n){
            rho=LANORM(&v[initer],2); 
            Jgrot[0][initer]=v[initer]/rho;
            Jgrot[1][initer]=v[initer+1]/rho;
            wgrot[initer+1]=-Jgrot[1][initer]*wgrot[initer];
            wgrot[initer]=conj(Jgrot[0][initer])*wgrot[initer];
            v[initer]=rho;
            v[initer+1]=complex<double>(0.0, 0.0);
        }

        for(int i=0; i<=initer; i++)
            Rgrot[i][initer]=v[i];
        normr=abs(wgrot[initer+1]);
        resvec[initer+1]=normr;
        normr_act=normr;
        if (pc_in==0) printf("GMRES-Iter.: %d, Res. Error: %1.3e\n", initer,normr/resvec[0]);std::cout.flush();
        

        if((normr<=tolb)||(initer==(inner-1))){

            complex<double>* ytmp=SOLVELSP(Rgrot,wgrot,initer);
            for(int i=0;i<n;i++)
                additive[i]=Ugrot[initer][i]*(-2.0*ytmp[initer]*conj(Ugrot[initer][initer]));
            additive[initer]+=ytmp[initer];
    
            for(int k=(initer-1); k>=0; k--){
                additive[k]+=ytmp[k];
                vec_temp=&Ugrot[k][0];
                val_product=2.0*SCALARP(vec_temp,additive,n);
                for(int i=0; i<n; i++){
                    additive[i]-=Ugrot[k][i]*val_product;
                }
            }
            for(int i=0; i<n; i++){
                vec_sol[i]+=additive[i];
            }
            

            if (pc_in!=0) printf("GMRES-IN-Iter.: %d, Res. Error: %1.3e\n", initer,normr/resvec[0]);

            destructor2D(Jgrot, 2, inner);
            destructor2D(Ugrot, inner+1, n);
            destructor2D(Rgrot, inner, inner);
            delete[] wgrot;
            delete[] Uhouse;
            delete[] v;
            delete[] vec_mvm;
            delete[] r;
            delete[] additive;
            delete[] ytmp;

            return 0;
        }
        
    }
    

    std::cout.flush();
    return 1;
}


double LANORM(complex<double> vec_in[],int size_vec){
    // 2-norm of a 1-D complex array "vec_in" with size "size_vec"
    double val_norm = 0.0;
    complex<double> val_in;
    for(int i=0; i<size_vec; i++){
        val_in=vec_in[i];
        val_norm += val_in.real()*val_in.real() + val_in.imag()*val_in.imag();
    }

    val_norm=sqrt(val_norm);
    return val_norm;
}

void cleartozero(complex<double> vec_in[],int size_vec){
    // Set the values of a 1-D complex array "vec_in" with size "size_vec" to zero
    for(int i=0; i<size_vec; i++)
        vec_in[i]=complex<double>(0.0, 0.0);

}

complex<double> ** constructor2D(int num_row, int num_col){
    complex<double> **data_in = new complex<double>*[num_row];
    for(int i=0; i<num_row; i++){
        data_in[i]=new complex<double>[num_col];
    }

    for(int i=0; i<num_row; i++){
        for(int j=0; j<num_col; j++)
            data_in[i][j]=complex<double>(0.0, 0.0);
    }

    return data_in;
}

void destructor2D(complex<double> **data_in, int num_row, int num_col){
    for(int i=0; i<num_row; i++){
        delete[] data_in[i];
    }
    delete[] data_in;
}

complex<double> SCALARSIGN(complex<double> inputval){
    complex<double> signofinput;
    if(abs(inputval)==0.0){
        signofinput=complex<double>(1.0, 0.0);
    } else {
        signofinput=inputval/abs(inputval);
    }
    return signofinput;
}

complex<double> SCALARP(complex<double>* vec_in1,complex<double>* vec_in2,int size_vec){
    complex<double> multiplication=complex<double>(0.0, 0.0);
    for(int i=0; i<size_vec; i++){
        multiplication+=conj(vec_in1[i])*vec_in2[i];
    }
    return multiplication;
}


complex<double>* SOLVELSP(complex<double>** A, complex<double>* b, int sizeA){
    complex<double>* x=new complex<double>[sizeA+1];
    for(int i=(sizeA); i>=0; i--){
        for(int j=(sizeA); j>=(i+1); j--){
            b[i]-=A[i][j]*x[j];
        }
        x[i]=b[i]/A[i][i];
    }
    return x;
}