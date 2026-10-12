#include <stdio.h>
#include <complex>
#include <chrono>
#include <iostream>
//#include <omp.h>

using namespace std;
using namespace std::chrono;

extern double pi;
extern double val_etaO;
extern double val_cfiealpha;
extern double val_kO;
extern int source_type;
extern complex<double> *source_CSP;
extern complex<double> *source_current;
extern complex<double> *source_DG;
extern double *source_dipole;
extern double *source_PW;
extern double val_epsilonO;
extern double *vec_trinormals;
extern double val_muO;
extern double *node_pos;
extern int num_rwgs;
extern int num_tris;
extern double val_omega;
extern complex<double> *vec_RHS;
extern double *rwg_lengths;
extern double *tri_areas;
extern int *tri_nodes;
extern int *tri2rwgindex;
extern int *tri2rwgnode;
extern int *tri2rwgsign;

double dotProduct3d(double vec1[], double vec2[]);
complex<double> dotProduct3d_dc(double vec1[], complex<double> vec2[]);
void crossProduct3d(double vec1[], double vec2[], double resvec[]);
void crossProduct3d_dc(double vec1[], complex<double> vec2[], complex<double> resvec[]);


void calculate_rhsvector(double dir_PW[]){
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int testnode1,testnode2,testnode3;
    int testfunction;
    int testnode;
    double testsign,testlength;
    double testnode1loc[3],testnode2loc[3],testnode3loc[3];
    double testlocation[3];
    double trinormal[3];
    double dotval;
    double vec_temp1[3];
    complex<double> vec_cross[3];
    complex<double> efield[3];
    complex<double> hfield[3],nhfield[3];
    complex<double> integral_efield,integral_hfield;
    complex<double> i=complex<double>(0.0, 1.0);
    complex<double> complex0=complex<double>(0.0, 0.0);

    vec_RHS=new complex<double>[num_rwgs*2];
    for(int ind0=0;ind0<(num_rwgs*2);ind0++){
        vec_RHS[ind0]=complex0;
    }

    printf("RHS Vector Size: %d  [%f GB]\n", num_rwgs, num_rwgs*16.0/1024.0/1024/1024);

    //RHS vector for plane wave illumination
    if(source_type==0){
         double pw_rad[4]={dir_PW[0]*pi/180.0,dir_PW[1]*pi/180.0,dir_PW[2]*pi/180.0,dir_PW[3]*pi/180.0};
         double pw_cart_pro[3]={sin(pw_rad[0])*cos(pw_rad[1]), sin(pw_rad[0])*sin(pw_rad[1]), cos(pw_rad[0])};
         double pw_cart_pol[3]={sin(pw_rad[2])*cos(pw_rad[3]), sin(pw_rad[2])*sin(pw_rad[3]), cos(pw_rad[2])};

         printf("Illumination: PW(deg) %1.2f %1.2f %1.2f %1.2f \n", dir_PW[0],dir_PW[1],dir_PW[2],dir_PW[3]);

         printf("Illumination: pw_cart_pro %1.2f %1.2f %1.2f  \n", pw_cart_pro[0],pw_cart_pro[1],pw_cart_pro[2]);
         printf("Illumination: pw_cart_pol %1.2f %1.2f %1.2f  \n", pw_cart_pol[0],pw_cart_pol[1],pw_cart_pol[2]);

        for(int testtri=0;testtri<num_tris;testtri++){
            testnode1=tri_nodes[testtri*3+0];
            testnode2=tri_nodes[testtri*3+1];
            testnode3=tri_nodes[testtri*3+2];
            testnode1loc[0]=node_pos[testnode1*3+0];
            testnode1loc[1]=node_pos[testnode1*3+1];
            testnode1loc[2]=node_pos[testnode1*3+2];
            testnode2loc[0]=node_pos[testnode2*3+0];
            testnode2loc[1]=node_pos[testnode2*3+1];
            testnode2loc[2]=node_pos[testnode2*3+2];
            testnode3loc[0]=node_pos[testnode3*3+0];
            testnode3loc[1]=node_pos[testnode3*3+1];
            testnode3loc[2]=node_pos[testnode3*3+2];
            testlocation[0]=(testnode1loc[0]+testnode2loc[0]+testnode3loc[0])/3;
            testlocation[1]=(testnode1loc[1]+testnode2loc[1]+testnode3loc[1])/3;
            testlocation[2]=(testnode1loc[2]+testnode2loc[2]+testnode3loc[2])/3;
            trinormal[0]=vec_trinormals[testtri*3+0];
            trinormal[1]=vec_trinormals[testtri*3+1];
            trinormal[2]=vec_trinormals[testtri*3+2];
            dotval=dotProduct3d(pw_cart_pro,testlocation);
            efield[0]=pw_cart_pol[0]*exp(i*val_kO*dotval);
            efield[1]=pw_cart_pol[1]*exp(i*val_kO*dotval);
            efield[2]=pw_cart_pol[2]*exp(i*val_kO*dotval);
            crossProduct3d_dc(pw_cart_pro,efield,vec_cross);
            hfield[0]=vec_cross[0]/val_etaO;
            hfield[1]=vec_cross[1]/val_etaO;
            hfield[2]=vec_cross[2]/val_etaO;
            crossProduct3d_dc(trinormal,hfield,nhfield);
            for(int indedge=0;indedge<3;indedge++){
                testfunction=tri2rwgindex[testtri*3+indedge];
                if(testfunction!=-1){
                    testnode=tri2rwgnode[testtri*3+indedge];
                    testsign=double(tri2rwgsign[testtri*3+indedge]);
                    testlength=rwg_lengths[testfunction];
                    vec_temp1[0]=testlocation[0]-node_pos[testnode*3+0];
                    vec_temp1[1]=testlocation[1]-node_pos[testnode*3+1];
                    vec_temp1[2]=testlocation[2]-node_pos[testnode*3+2];
                    integral_efield=dotProduct3d_dc(vec_temp1,efield);
                    integral_hfield=dotProduct3d_dc(vec_temp1,hfield);
                    vec_RHS[testfunction]-=testsign*testlength*integral_efield/2.0;
                    vec_RHS[testfunction+num_rwgs]-=testsign*testlength*integral_hfield/2.0*val_etaO;
                }
            }
        }
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    printf("RHS Vector Calculation Takes %f s \n", time_span.count());

    std::cout.flush();
}

