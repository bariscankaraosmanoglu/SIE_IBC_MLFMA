#include <stdio.h>
#include <cmath>
#include <complex>
#include <chrono>
#include <iostream>

using namespace std;

extern int num_tris;
extern int num_rwgs;
extern int *rwg_nodes;
extern int *rwg_tris;
extern double *node_pos;
extern double *rwg_lengths;
extern complex<double> *vec_sol;
extern double dir_cfc[];
extern double val_kO;
extern double val_muO;
extern double val_omega;
extern double val_etaO;
extern int *num_filledboxes;
extern int val_maxtri;
extern int *smallestbox2tri;
extern int num_mlfmalevels;
extern double pi;


void calculate_currentfieldcontribution(){

    int indtri;
    int rwglastnodeindex[2]{2,3};
    int rwgnode1,rwgnode2,rwgnode3;
    double rwgsignval;
    double rwgsign[2]{1.0,-1.0};
    double rwgnode1loc[3],rwgnode2loc[3],rwgnode3loc[3],rwglocation[3];
    double rwglength;
    double khat[3];
    double val_dot;
    complex<double> rwg_coef,rwgvector[3];
    complex<double> *farfieldoftris;
    double *cfcontris;
    complex<double> complex0=complex<double>(0.0, 0.0);
    complex<double> i=complex<double>(0.0, 1.0);
    complex<double> farfieldtotal[3]{complex0,complex0,complex0};
    double farfieldtotalnorm;

    farfieldoftris=new complex<double>[num_tris*3];
    for(int ind=0;ind<(num_tris*3);ind++){
        farfieldoftris[ind]=complex0;
    }
    cfcontris=new double[num_tris];
    for(int ind=0;ind<num_tris;ind++){
        cfcontris[ind]=0.0;
    }

    khat[0]=sin(dir_cfc[0]*pi/180.0)*cos(dir_cfc[1]*pi/180.0);
    khat[1]=sin(dir_cfc[0]*pi/180.0)*sin(dir_cfc[1]*pi/180.0);
    khat[2]=cos(dir_cfc[0]*pi/180.0);

    for(int indrwg=0;indrwg<num_rwgs;indrwg++){
        for(int indrwgtri=0;indrwgtri<2;indrwgtri++){
            indtri=rwg_tris[indrwg*2+indrwgtri];
            rwgsignval=rwgsign[indrwgtri];
            rwgnode1=rwg_nodes[indrwg*4+0];
            rwgnode2=rwg_nodes[indrwg*4+1];
            rwgnode3=rwg_nodes[indrwg*4+rwglastnodeindex[indrwgtri]];
            rwgnode1loc[0]=node_pos[rwgnode1*3+0];
            rwgnode1loc[1]=node_pos[rwgnode1*3+1];
            rwgnode1loc[2]=node_pos[rwgnode1*3+2];
            rwgnode2loc[0]=node_pos[rwgnode2*3+0];
            rwgnode2loc[1]=node_pos[rwgnode2*3+1];
            rwgnode2loc[2]=node_pos[rwgnode2*3+2];
            rwgnode3loc[0]=node_pos[rwgnode3*3+0];
            rwgnode3loc[1]=node_pos[rwgnode3*3+1];
            rwgnode3loc[2]=node_pos[rwgnode3*3+2];
            rwglocation[0]=(rwgnode1loc[0]+rwgnode2loc[0]+rwgnode3loc[0])/3;
            rwglocation[1]=(rwgnode1loc[1]+rwgnode2loc[1]+rwgnode3loc[1])/3;
            rwglocation[2]=(rwgnode1loc[2]+rwgnode2loc[2]+rwgnode3loc[2])/3;
            val_dot=khat[0]*rwglocation[0]+khat[1]*rwglocation[1]+khat[2]*rwglocation[2];
            rwglength=rwg_lengths[indrwg];
            rwgvector[0]=rwgsignval*rwglength*(rwglocation[0]-rwgnode3loc[0]); //single point integration
            rwgvector[1]=rwgsignval*rwglength*(rwglocation[1]-rwgnode3loc[1]);
            rwgvector[2]=rwgsignval*rwglength*(rwglocation[2]-rwgnode3loc[2]);
            rwg_coef=vec_sol[indrwg]/val_etaO;
            farfieldoftris[indtri*3+0]+=i*val_omega*val_muO*exp(-i*val_kO*val_dot)/(8.0*pi)*rwgvector[0]*rwg_coef;
            farfieldoftris[indtri*3+1]+=i*val_omega*val_muO*exp(-i*val_kO*val_dot)/(8.0*pi)*rwgvector[1]*rwg_coef;
            farfieldoftris[indtri*3+2]+=i*val_omega*val_muO*exp(-i*val_kO*val_dot)/(8.0*pi)*rwgvector[2]*rwg_coef;

            complex<double> temp_farfield[3];
            rwg_coef=vec_sol[indrwg+num_rwgs];
            temp_farfield[0]=-i*val_kO*exp(-i*val_kO*val_dot)/(8.0*pi)*rwgvector[0]*rwg_coef;
            temp_farfield[1]=-i*val_kO*exp(-i*val_kO*val_dot)/(8.0*pi)*rwgvector[1]*rwg_coef;
            temp_farfield[2]=-i*val_kO*exp(-i*val_kO*val_dot)/(8.0*pi)*rwgvector[2]*rwg_coef;
            farfieldoftris[indtri*3+0]+=khat[1]*temp_farfield[2]-khat[2]*temp_farfield[1];
            farfieldoftris[indtri*3+1]+=khat[2]*temp_farfield[0]-khat[0]*temp_farfield[2];
            farfieldoftris[indtri*3+2]+=khat[0]*temp_farfield[1]-khat[1]*temp_farfield[0];
        }
    }

    for(int ind=0;ind<num_tris;ind++){
        farfieldtotal[0]+=farfieldoftris[ind*3+0];
        farfieldtotal[1]+=farfieldoftris[ind*3+1];
        farfieldtotal[2]+=farfieldoftris[ind*3+2];
    }
    farfieldtotalnorm=sqrt(abs(farfieldtotal[0])*abs(farfieldtotal[0])+\
                           abs(farfieldtotal[1])*abs(farfieldtotal[1])+\
                           abs(farfieldtotal[2])*abs(farfieldtotal[2]));

    printf("farfieldtotalnorm: %1.8e \n", farfieldtotalnorm);


    int box_start=num_filledboxes[num_mlfmalevels*2-1];
    int box_end=num_filledboxes[num_mlfmalevels*3-1]+1;
    int basistri;
    double farfieldofboxrratio;
    for(int box_basis=box_start;box_basis<box_end;box_basis++){
        complex<double> farfieldofbox[3]{complex0,complex0,complex0};
        for(int indbasistri=0;indbasistri<val_maxtri;indbasistri++){
            basistri=smallestbox2tri[(box_basis-box_start)*val_maxtri+indbasistri];
            if(basistri!=-1){
                farfieldofbox[0]+=farfieldoftris[basistri*3+0];
                farfieldofbox[1]+=farfieldoftris[basistri*3+1];
                farfieldofbox[2]+=farfieldoftris[basistri*3+2];
            }
        }
        farfieldofboxrratio=sqrt(abs(farfieldofbox[0])*abs(farfieldofbox[0])+\
                                 abs(farfieldofbox[1])*abs(farfieldofbox[1])+\
                                 abs(farfieldofbox[2])*abs(farfieldofbox[2]))/farfieldtotalnorm;
        for(int indbasistri=0;indbasistri<val_maxtri;indbasistri++){
            basistri=smallestbox2tri[(box_basis-box_start)*val_maxtri+indbasistri];
            if(basistri!=-1){
                cfcontris[basistri]=farfieldofboxrratio;
            }
        }
    }


    FILE *fptr;
    fptr=fopen("result_cfcontris.txt","w");
    for(int ind=0;ind<num_tris;ind++){
        fprintf(fptr,"%1.8e\n", cfcontris[ind]);
    }
    fclose(fptr);


    delete[] farfieldoftris;
    delete[] cfcontris;

    std::cout.flush();

}
