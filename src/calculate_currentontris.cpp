#include <stdio.h>
#include <cmath>
#include <complex>
#include <chrono>
#include <iostream>

using namespace std;
using namespace std::chrono;


extern int num_tris;
extern int num_rwgs;
extern int *rwg_nodes;
extern int *rwg_tris;
extern double *node_pos;
extern double *rwg_lengths;
extern double *tri_areas;
extern double val_etaO;
extern complex<double> *vec_sol;

void calculate_currentontris(){

    int indtri;
    int rwglastnodeindex[2]{2,3};
    int rwgnode1,rwgnode2,rwgnode3;
    double rwgsignval;
    double rwgsign[2]{1.0,-1.0};
    double rwgnode1loc[3],rwgnode2loc[3],rwgnode3loc[3],rwglocation[3];
    double rwglength;
    double tri_area;
    complex<double> rwg_coef,rwgvector[3];
    complex<double> *currentontris;
    complex<double> complex0=complex<double>(0.0, 0.0);

    currentontris=new complex<double>[num_tris*3];
    for(int ind=0;ind<(num_tris*3);ind++){
        currentontris[ind]=complex0;
    }

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
            rwg_coef=vec_sol[indrwg]/val_etaO;
            rwglength=rwg_lengths[indrwg];
            tri_area=tri_areas[indtri];
            rwgvector[0]=rwgsignval*rwg_coef*rwglength/tri_area*(rwglocation[0]-rwgnode3loc[0]);
            rwgvector[1]=rwgsignval*rwg_coef*rwglength/tri_area*(rwglocation[1]-rwgnode3loc[1]);
            rwgvector[2]=rwgsignval*rwg_coef*rwglength/tri_area*(rwglocation[2]-rwgnode3loc[2]);
            currentontris[indtri*3+0]+=rwgvector[0];
            currentontris[indtri*3+1]+=rwgvector[1];
            currentontris[indtri*3+2]+=rwgvector[2];
        }
    }


    FILE *fptr;
    fptr=fopen("result_currentontris_J.txt","w");
    for(int ind=0;ind<num_tris;ind++){
        fprintf(fptr,"%1.8e\n", sqrt(abs(currentontris[ind*3+0])*abs(currentontris[ind*3+0])+ \
                                     abs(currentontris[ind*3+1])*abs(currentontris[ind*3+1])+ \
                                     abs(currentontris[ind*3+2])*abs(currentontris[ind*3+2])));
    }
    fclose(fptr);



    for(int ind=0;ind<(num_tris*3);ind++){
        currentontris[ind]=complex0;
    }

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
            rwg_coef=vec_sol[indrwg+num_rwgs];
            rwglength=rwg_lengths[indrwg];
            tri_area=tri_areas[indtri];
            rwgvector[0]=rwgsignval*rwg_coef*rwglength/tri_area*(rwglocation[0]-rwgnode3loc[0]);
            rwgvector[1]=rwgsignval*rwg_coef*rwglength/tri_area*(rwglocation[1]-rwgnode3loc[1]);
            rwgvector[2]=rwgsignval*rwg_coef*rwglength/tri_area*(rwglocation[2]-rwgnode3loc[2]);
            currentontris[indtri*3+0]+=rwgvector[0];
            currentontris[indtri*3+1]+=rwgvector[1];
            currentontris[indtri*3+2]+=rwgvector[2];
        }
    }


    fptr=fopen("result_currentontris_M.txt","w");
    for(int ind=0;ind<num_tris;ind++){
        fprintf(fptr,"%1.8e\n", sqrt(abs(currentontris[ind*3+0])*abs(currentontris[ind*3+0])+ \
                                     abs(currentontris[ind*3+1])*abs(currentontris[ind*3+1])+ \
                                     abs(currentontris[ind*3+2])*abs(currentontris[ind*3+2])));
    }
    fclose(fptr);

    delete[] currentontris;


    std::cout.flush();


}
