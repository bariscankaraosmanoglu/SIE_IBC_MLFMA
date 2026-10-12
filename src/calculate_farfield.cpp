#include "calculate_farfield.h"

void calculate_farfield(int ind_sol){
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    //
    //
    // calculate_farfield
    // Calculating the scattered electric fields from equivalent currents
    double loc_theta=0;
    double loc_phi=0;
    double costheta=0;
    double sintheta=0;
    double cosphi=0;
    double sinphi=0;
    double khat[3]{0,0,0};
    int basisnode;
    int basissign;
    int indbasis;
    int basisnodes[3];
    double basislength;
    double basislocation[3]{0,0,0};
    double basisnode0loc[3]{0,0,0};
    double basisnode1loc[3]{0,0,0};
    double basisnode2loc[3]{0,0,0};
    double basisnode3loc[3]{0,0,0};
    int num_quadpnts=3;
    double basisweight, basisweights[num_quadpnts], basislocations[num_quadpnts][3];
    double val_dot;
    complex<double> val_i=complex<double>(0.0, 1.0);
    complex<double> temp_farfield[3];
    complex<double> temp_farfield_sum[3];
    complex<double> val_product1=complex<double>(0.0, 0.0);
    complex<double> val_product2=complex<double>(0.0, 0.0);
    complex<double> val_product3=complex<double>(0.0, 0.0);
    complex<double> val_theta;
    complex<double> val_phi;
    FILE *fptr;


    if(dir_ff==0){
        num_ffdir=1;
        ffdir_theta = new double[num_ffdir];
        ffdir_phi = new double[num_ffdir];
        ffdir_theta[0]=fmod(180-dir_PW[0], 180)*pi/180.0;
        if(dir_PW[0]==0) ffdir_theta[0]=pi;
        ffdir_phi[0]=fmod(dir_PW[1]+180, 360)*pi/180.0;
    }
    else if(dir_ff>0){
        if (dir_ff<0.1) dir_ff=0.1;
        double del_angle=dir_ff;
        num_ffdir=(180/del_angle + 1)*(360/del_angle + 1);
        ffdir_theta = new double[num_ffdir];
        ffdir_phi = new double[num_ffdir];  
        int ind_ff=0;
        for(int ind_theta=0; ind_theta<(180/del_angle+1); ind_theta++){
            for(int ind_phi=0; ind_phi<(360/del_angle+1); ind_phi++){
                ffdir_theta[ind_ff]=(del_angle*ind_theta)*pi/180.0;
                ffdir_phi[ind_ff]=(del_angle*ind_phi)*pi/180.0;
                ind_ff++;
            }
        }
    }
    else{
        num_ffdir_cut = 361;
        num_ffdir = num_ffdir_cut*3;
        ffdir_phi = new double[num_ffdir];
        ffdir_theta = new double[num_ffdir];
        
        for(int ind1=0;ind1<num_ffdir;ind1++){
            ffdir_phi[ind1]=0;
            ffdir_theta[ind1]=0;
        }
        int ind0=0;
        // Setting theta and phi angles on ZX-plane
        for(double ind1=0;ind1<(num_ffdir_cut-1)/2;ind1++){
            ffdir_theta[ind0]=pi*ind1/(num_ffdir_cut-1)*2;
            ffdir_phi[ind0]=0;
            ind0++;
        }
        for(double ind1=(num_ffdir_cut-1)/2;ind1<num_ffdir_cut;ind1++){
            ffdir_theta[ind0]=2*pi-pi*(ind1/(num_ffdir_cut-1))*2;
            ffdir_phi[ind0]=pi;
            ind0++;
        }
        // Setting theta and phi angles on ZY-plane
        for(double ind1=0;ind1<(num_ffdir_cut-1)/2;ind1++){
            ffdir_theta[ind0]=pi*ind1/(num_ffdir_cut-1)*2;
            ffdir_phi[ind0]=pi/2;
            ind0++;
        }
        for(double ind1=(num_ffdir_cut-1)/2;ind1<num_ffdir_cut;ind1++){
            ffdir_theta[ind0]=2*pi-pi*(ind1/(num_ffdir_cut-1))*2;
            ffdir_phi[ind0]=3*pi/2;
            ind0++;
        }
        // Setting theta and phi angles on XY-plane
        for(double ind1=0;ind1<num_ffdir_cut;ind1++){
            ffdir_theta[ind0]=pi/2;
            ffdir_phi[ind0]=pi*ind1/(num_ffdir_cut-1)*2;
            ind0++;
        }
    }

    if(dir_ff==0){
        if(ind_sol==0){
            fptr=fopen("result_farefield.txt","w");
            fprintf(fptr,"%%Theta       Phi         Abs(E_theta)    Phase(E_theta)  Abs(E_phi)      Phase(E_phi)\n");
        }
        else{
            fptr=fopen("result_farefield.txt","a");
        }
    } else{
        char filename[100];
        sprintf(filename, "result_farefield_%d.txt", ind_sol);
        fptr=fopen(filename,"w");
        fprintf(fptr,"%%Theta       Phi         Abs(E_theta)    Phase(E_theta)  Abs(E_phi)      Phase(E_phi)\n");
    }



    vec_farfield=new complex<double>[num_ffdir*2];
    for(int ind1=0;ind1<(num_ffdir*2);ind1++)
        vec_farfield[ind1]=complex<double>(0.0, 0.0);
    for(int indloc=0; indloc<num_ffdir; indloc++){
        temp_farfield_sum[0]=complex<double>(0.0, 0.0);
        temp_farfield_sum[1]=complex<double>(0.0, 0.0);
        temp_farfield_sum[2]=complex<double>(0.0, 0.0);
        loc_theta=ffdir_theta[indloc];
        loc_phi=ffdir_phi[indloc];
        costheta=cos(loc_theta);
        sintheta=sin(loc_theta);
        cosphi=cos(loc_phi);
        sinphi=sin(loc_phi);
        khat[0]=sintheta*cosphi;
        khat[1]=sintheta*sinphi;
        khat[2]=costheta;
        for(int basistri=0; basistri<num_tris; basistri++){
            basisnodes[0]=tri_nodes[basistri*3+0];
            basisnodes[1]=tri_nodes[basistri*3+1];
            basisnodes[2]=tri_nodes[basistri*3+2];
            basisnode1loc[0]=node_pos[basisnodes[0]*3+0];
            basisnode1loc[1]=node_pos[basisnodes[0]*3+1];
            basisnode1loc[2]=node_pos[basisnodes[0]*3+2];
            basisnode2loc[0]=node_pos[basisnodes[1]*3+0];
            basisnode2loc[1]=node_pos[basisnodes[1]*3+1];
            basisnode2loc[2]=node_pos[basisnodes[1]*3+2];
            basisnode3loc[0]=node_pos[basisnodes[2]*3+0];
            basisnode3loc[1]=node_pos[basisnodes[2]*3+1];
            basisnode3loc[2]=node_pos[basisnodes[2]*3+2];
            quadrature_points(basisnode1loc,basisnode2loc,basisnode3loc,basisweights,basislocations,num_quadpnts);
            for(int ind1=0; ind1<3; ind1++){
                indbasis=tri2rwgindex[basistri*3+ind1];
                if(indbasis!=-1){
                    basisnode=tri2rwgnode[basistri*3+ind1];
                    basissign=tri2rwgsign[basistri*3+ind1];
                    basislength=rwg_lengths[indbasis];
                    basisnode0loc[0]=node_pos[basisnode*3+0];
                    basisnode0loc[1]=node_pos[basisnode*3+1];
                    basisnode0loc[2]=node_pos[basisnode*3+2];
                    for(int indb=0;indb<num_quadpnts;indb++){
                        basisweight=basisweights[indb];
                        basislocation[0]=basislocations[indb][0];
                        basislocation[1]=basislocations[indb][1];
                        basislocation[2]=basislocations[indb][2];
                        val_dot=khat[0]*basislocation[0]+khat[1]*basislocation[1]+khat[2]*basislocation[2];
                        val_product1=exp(-val_i*(val_kO*val_dot));
                        temp_farfield[0]=(basislocation[0]-basisnode0loc[0])*val_product1;
                        temp_farfield[1]=(basislocation[1]-basisnode0loc[1])*val_product1;
                        temp_farfield[2]=(basislocation[2]-basisnode0loc[2])*val_product1;
                        val_product2=val_i*val_omega*val_muO*basislength/(8*pi)*(vec_sol[indbasis]/val_etaO)*basisweight;
                        val_product2*=basissign;
                        val_product3=-val_i*val_kO*basislength/(8*pi)*vec_sol[indbasis+num_rwgs]*basisweight;
                        val_product3*=basissign;
                        complex<double> temp_farfield_cross[3];
                        temp_farfield_cross[0]=khat[1]*temp_farfield[2]-khat[2]*temp_farfield[1];
                        temp_farfield_cross[1]=khat[2]*temp_farfield[0]-khat[0]*temp_farfield[2];
                        temp_farfield_cross[2]=khat[0]*temp_farfield[1]-khat[1]*temp_farfield[0];

                        temp_farfield_sum[0]+=temp_farfield[0]*val_product2;
                        temp_farfield_sum[1]+=temp_farfield[1]*val_product2;
                        temp_farfield_sum[2]+=temp_farfield[2]*val_product2;

                        temp_farfield_sum[0]+=temp_farfield_cross[0]*val_product3;
                        temp_farfield_sum[1]+=temp_farfield_cross[1]*val_product3;
                        temp_farfield_sum[2]+=temp_farfield_cross[2]*val_product3;

                    }

                }
            }
        }
        val_theta=costheta*cosphi*temp_farfield_sum[0]+costheta*sinphi*temp_farfield_sum[1]-sintheta*temp_farfield_sum[2];
        val_phi=-sinphi*temp_farfield_sum[0]+cosphi*temp_farfield_sum[1];
        vec_farfield[indloc*2+0]=val_theta;
        vec_farfield[indloc*2+1]=val_phi;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    printf("Farfield Calculation Takes %f s \n", time_span.count());


    //fptr=fopen("result_farefield.txt","w");
    //fprintf(fptr,"%%Theta       Phi         Abs(E_theta)    Phase(E_theta)  Abs(E_phi)      Phase(E_phi)\n");
    for(int ind=0;ind<num_ffdir;ind++){
        fprintf(fptr,"%1.4e %1.4e %1.8e %1.8e %1.8e %1.8e \n", ffdir_theta[ind]*180.0/pi, ffdir_phi[ind]*180.0/pi, \
                 abs(vec_farfield[ind*2+0]), arg(vec_farfield[ind*2+0]), \
                 abs(vec_farfield[ind*2+1]), arg(vec_farfield[ind*2+1]));
    }
    fclose(fptr);
    delete[] ffdir_theta; ffdir_theta = nullptr;
    delete[] ffdir_phi; ffdir_phi = nullptr;
    delete[] vec_farfield; vec_farfield = nullptr;


    std::cout.flush();

}
