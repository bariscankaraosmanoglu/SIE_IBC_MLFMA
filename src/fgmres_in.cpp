#include <stdio.h>
#include <complex>
#include <cmath>
#include <iostream>

using namespace std;

void mvm_in(complex<double> x[], complex<double> y[]);
void pc_solution_in(complex<double> vec_RHS_i[], complex<double> vec_sol_i[], int size_vecRHS);

int fgmres_in(int size_vecRHS,complex<double> vec_RHS[],int itermax,int iterrestart,double soltolerance,complex<double> vec_sol[]){


    int its=0;
    int i=0;
    int i1=0;
    double normrhs=0.0;
    double ro=0.0;
    double tol1=0.0;
    double *res;


    complex<double> complex0=complex<double>(0.0, 0.0);
    complex<double> *vv;
    complex<double> *hh;
    complex<double> *Z;
    complex<double> *vec_mvm;
    complex<double> *z;
    complex<double> *vec_temp1;
    complex<double> t=complex<double>(0.0, 0.0);
    complex<double> *rs;
    complex<double> *c;
    complex<double> *s;
    complex<double> gam=complex0;


    for(int ind0=0;ind0<size_vecRHS;ind0++){
        normrhs+=vec_RHS[ind0].real()*vec_RHS[ind0].real()+vec_RHS[ind0].imag()*vec_RHS[ind0].imag();
    }
    normrhs=sqrt(normrhs);
//    printf("normrhs: %e \n", normrhs);
    vv=new complex<double>[size_vecRHS*(iterrestart+1)];
    hh=new complex<double>[(iterrestart+1)*(iterrestart+1)];
    Z=new complex<double>[size_vecRHS*(itermax+1)];
    res=new double[itermax+2];
    vec_mvm=new complex<double>[size_vecRHS];
    vec_temp1=new complex<double>[size_vecRHS];
    z=new complex<double>[size_vecRHS];
    rs=new complex<double>[iterrestart+1];
    c=new complex<double>[iterrestart+1];
    s=new complex<double>[iterrestart+1];

//    printf("FGMRES Solver Memory: %d  [%f GB]\n", size_vecRHS*(iterrestart+itermax+6), size_vecRHS*(iterrestart+itermax+6)*16.0/1024.0/1024/1024);

    for(int ind=0;ind<(size_vecRHS*iterrestart);ind++){
            vv[ind]=complex0;
    }
    for(int ind=0;ind<(iterrestart*iterrestart);ind++){
        hh[ind]=complex0;
    }
    for(int ind=0;ind<(size_vecRHS*itermax);ind++){
        Z[ind]=complex0;
    }
    for(int ind=0;ind<(size_vecRHS);ind++){
        vec_mvm[ind]=complex0;
        z[ind]=complex0;
    }
    for(int ind=0;ind<=iterrestart;ind++){
        rs[ind]=complex0;
        c[ind]=complex0;
        s[ind]=complex0;
    }

    while(its<itermax){
        mvm_in(vec_sol,vec_mvm);
        for(int ind=0;ind<size_vecRHS;ind++){
            vv[ind]=vec_RHS[ind]-vec_mvm[ind];
        }
        double val_norm=0.0;
        for(int ind=0;ind<size_vecRHS;ind++){
            val_norm+=vv[ind].real()*vv[ind].real()+vv[ind].imag()*vv[ind].imag();
        }
        val_norm=sqrt(val_norm);
        ro=val_norm;
        res[its]=ro/normrhs;
        if(its==0){
            tol1=soltolerance*normrhs;
        }
        if((ro<=tol1)||(its>=itermax)){
            return 0;
        }
        t=1.0/ro;
        for(int ind=0;ind<size_vecRHS;ind++){
            vv[ind]*=t;
        }
        rs[0]=ro;
//        printf("FGMRES-IN-Iter.: %d, Res. Error: %1.3e \n", its,ro/normrhs);
        i=0;

        while((i<iterrestart)&&(ro>tol1)&&(its<itermax)){
            i1=i+1;
            for(int ind=0;ind<size_vecRHS;ind++){
                vec_temp1[ind]=vv[i*size_vecRHS+ind];
            }
            pc_solution_in(vec_temp1,z,size_vecRHS);
            for(int ind=0;ind<size_vecRHS;ind++){
                Z[i*size_vecRHS+ind]=z[ind];
            }
            for(int ind=0;ind<(size_vecRHS);ind++){
                vec_mvm[ind]=complex0;
            }
            mvm_in(z,vec_mvm);
            for(int ind=0;ind<size_vecRHS;ind++){
                vv[i1*size_vecRHS+ind]=vec_mvm[ind];
            }
            for(int j=0;j<=i;j++){
                t=complex0;
                for(int ind=0;ind<size_vecRHS;ind++){
                    t+=conj(vv[j*size_vecRHS+ind])*vv[i1*size_vecRHS+ind];
                }
                 hh[j*iterrestart+i]=t;
                for(int ind=0;ind<size_vecRHS;ind++){
                    vv[i1*size_vecRHS+ind]-=t*vv[j*size_vecRHS+ind];
                }
            }
            t=complex0;
            for(int ind=0;ind<size_vecRHS;ind++){
                t+=vv[i1*size_vecRHS+ind].real()*vv[i1*size_vecRHS+ind].real()+vv[i1*size_vecRHS+ind].imag()*vv[i1*size_vecRHS+ind].imag();
            }
            t=sqrt(t);
            hh[i1*iterrestart+i]=t;
            if(t!=complex0){
                t=1.0/t;
                for(int ind=0;ind<size_vecRHS;ind++){
                    vv[i1*size_vecRHS+ind]=t*vv[i1*size_vecRHS+ind];
                }
            }
            if(i!=0){
                for(int k=1;k<=i;k++){
                    int k1=k-1;
                    t=hh[k1*iterrestart+i];
                    hh[k1*iterrestart+i]=c[k1]*t+s[k1]*hh[k*iterrestart+i];
                    hh[k*iterrestart+i]=-s[k1]*t+c[k1]*hh[k*iterrestart+i];
                }
            }
            gam=sqrt(hh[i*iterrestart+i]*hh[i*iterrestart+i]+hh[i1*iterrestart+i]*hh[i1*iterrestart+i]);
            if(gam==complex0){
                gam=2.2e-16;
            }
            c[i]=hh[i*iterrestart+i]/gam;
            s[i]=hh[i1*iterrestart+i]/gam;
            rs[i1]=-s[i]*rs[i];
            rs[i]=c[i]*rs[i];
            hh[i*iterrestart+i]=c[i]*hh[i*iterrestart+i]+s[i]*hh[i1*iterrestart+i];
            ro=abs(rs[i1]);
            res[its+1]=ro;
//            printf("FGMRES-IN-Iter.: %d, Res. Error: %1.3e \n", its+1,ro/normrhs);
            i++;
            its++;
        }
        i--;

        rs[i]/=hh[(i)*iterrestart+(i)];
        for(int k=i-1;k>=0;k--){
            t=rs[k];
            for(int j=k+1;j<=i;j++){
                t-=hh[k*iterrestart+j]*rs[j];
            }
            rs[k]=t/hh[k*iterrestart+k];
        }
        for(int j=0;j<=i;j++){
            for(int ind=0;ind<size_vecRHS;ind++){
                vec_sol[ind]+=rs[j]*Z[j*size_vecRHS+ind];
            }
        }

        if((ro<=tol1) || (its>=itermax)){
            printf("FGMRES-IN-Iter.: %d, Res. Error: %1.3e \n", its+1,ro/normrhs);

            delete[] vv;
            delete[] hh;
            delete[] Z;
            delete[] vec_mvm;
            delete[] z;
            delete[] vec_temp1;
            delete[] rs;
            delete[] c;
            delete[] s;

            return 1;
        }
    }

    std::cout.flush();
    return 0;
}

