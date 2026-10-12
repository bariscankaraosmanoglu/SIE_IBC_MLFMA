#include "calculate_loadbalance.h"
#include "calculate_Znearmatrix.h"

void calculate_Znearmatrix(){
    //
    //
    // calculate_Znearmatrix
    // Calculate rwg-rwg interactions for the Znear matrix.
    // Triangle to triangle calculations are performed.

    int box_start=num_filledboxes[num_mlfmalevels*2-1];
    int box_end=num_filledboxes[num_mlfmalevels*3-1]+1;
    int testtri,basistri;
    int testnode1,testnode2,testnode3;
    int basisnode1,basisnode2,basisnode3;
    double testnode1loc[3],testnode2loc[3],testnode3loc[3];
    double basisnode1loc[3],basisnode2loc[3],basisnode3loc[3];
    double testtrinormal[3],basistrinormal[3];
    double testlocation[3];
    double basisarea;
    int num_quadpnts=3;
    double basisweight, basisweights[num_quadpnts], basislocations[num_quadpnts][3];
    double analytical1R[1],analyticalrhoR[3],analyticalgrad1R[3],analyticalgradR[3];
    double vec_R[3], distance;
    complex<double> i=complex<double>(0.0, 1.0);
    complex<double> complex0=complex<double>(0.0, 0.0);
    complex<double> extractedgreensAB, extractedgreensC;
    complex<double> integralA, integralB[3], integralC[3];
    int testfunction,testnode,testsign;
    double testlength;
    double vectorvalue2[3],test_node_pos[3],val_dot1,projectionvalue2[3];
    int basisfunction,basisnode,basissign;
    double basislength;
    double vectorvalue3[3],basis_node_pos[3],val_dot2,projectionvalue3[3];
    long int Znear_index;
    int box_basis;

    duration<double> time_span;
    high_resolution_clock::time_point t1,t2;

    t1 = high_resolution_clock::now();
    matrix_Znear11=new complex<double>[size_Znear];
    matrix_Znear12=new complex<double>[size_Znear];
    matrix_Znear21=new complex<double>[size_Znear];
    matrix_Znear22=new complex<double>[size_Znear];
    t2 = high_resolution_clock::now();
    time_span=t2-t1;
    //printf("Znear Allocation Takes %f s \n", time_span.count());


    t1 = high_resolution_clock::now();
    for(long int ind0=0;ind0<size_Znear;ind0++){
        matrix_Znear11[ind0]=complex0;
        matrix_Znear12[ind0]=complex0;
        matrix_Znear21[ind0]=complex0;
        matrix_Znear22[ind0]=complex0;
    }
    t2 = high_resolution_clock::now();
    time_span=t2-t1;
    //printf("Znear Initialization Takes %f s \n", time_span.count());
    printf("Znear Size: %ld  [%f GB]\n", size_Znear, size_Znear*16.0/1024.0/1024/1024);

    t1 = high_resolution_clock::now();
    omp_set_num_threads(num_threads);


    #pragma omp parallel private(box_basis,testtri,testnode1,testnode2,testnode3,testnode1loc,testnode2loc,testnode3loc,testtrinormal,testlocation,basistri,\
                                 basisnode1,basisnode2,basisnode3,basisnode1loc,basisnode2loc,basisnode3loc,basisarea,basistrinormal,basisweights,basislocations,\
                                 analytical1R,analyticalrhoR,analyticalgrad1R,analyticalgradR,integralA,integralB,integralC,basisweight,vec_R,distance,\
                                 extractedgreensAB,extractedgreensC,testfunction,testnode,testsign,testlength,vectorvalue2,test_node_pos,val_dot1,projectionvalue2,\
                                 basisfunction,basisnode,basissign,basislength,vectorvalue3,basis_node_pos,val_dot2,projectionvalue3,Znear_index) \
                                 if(parallel_omp_on)

                                 
    {
        #pragma omp barrier
        const auto box_range=near_box_range(omp_get_thread_num(),omp_get_num_threads());
        double timer_start = omp_get_wtime();
        double timer_total = 0.0;
        complex<double> val_imp{};

//        for(int box_test=box_start;box_test<box_end;box_test++){
        for(int box_test=box_range.first;box_test<box_range.second;box_test++){
            double timer_in = omp_get_wtime();
            for(int ind_nearbox=0;ind_nearbox<27;ind_nearbox++){
                if(boxes_near[box_test*27+ind_nearbox]!=-1){
                    box_basis=boxes_near[box_test*27+ind_nearbox];
                    for(int indtesttri=0;indtesttri<val_maxtri;indtesttri++){
                        testtri=smallestbox2tri[(box_test-box_start)*val_maxtri+indtesttri];
                        if(testtri!=-1){
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
                            testtrinormal[0]=vec_trinormals[testtri*3+0];
                            testtrinormal[1]=vec_trinormals[testtri*3+1];
                            testtrinormal[2]=vec_trinormals[testtri*3+2];
                            testlocation[0]=(testnode1loc[0]+testnode2loc[0]+testnode3loc[0])/3;
                            testlocation[1]=(testnode1loc[1]+testnode2loc[1]+testnode3loc[1])/3;
                            testlocation[2]=(testnode1loc[2]+testnode2loc[2]+testnode3loc[2])/3;
                            for(int indbasistri=0;indbasistri<val_maxtri;indbasistri++){
                                basistri=smallestbox2tri[(box_basis-box_start)*val_maxtri+indbasistri];
                                if(basistri!=-1){
                                    val_imp = val_impedance[tri_media[basistri]];
                                    basisnode1=tri_nodes[basistri*3+0];
                                    basisnode2=tri_nodes[basistri*3+1];
                                    basisnode3=tri_nodes[basistri*3+2];
                                    basisnode1loc[0]=node_pos[basisnode1*3+0];
                                    basisnode1loc[1]=node_pos[basisnode1*3+1];
                                    basisnode1loc[2]=node_pos[basisnode1*3+2];
                                    basisnode2loc[0]=node_pos[basisnode2*3+0];
                                    basisnode2loc[1]=node_pos[basisnode2*3+1];
                                    basisnode2loc[2]=node_pos[basisnode2*3+2];
                                    basisnode3loc[0]=node_pos[basisnode3*3+0];
                                    basisnode3loc[1]=node_pos[basisnode3*3+1];
                                    basisnode3loc[2]=node_pos[basisnode3*3+2];
                                    basisarea=tri_areas[basistri];
                                    basistrinormal[0]=vec_trinormals[basistri*3+0];
                                    basistrinormal[1]=vec_trinormals[basistri*3+1];
                                    basistrinormal[2]=vec_trinormals[basistri*3+2];
                                    quadrature_points(basisnode1loc,basisnode2loc,basisnode3loc,basisweights,basislocations,num_quadpnts);
                                    tri_analytical_int(basisnode1loc,basisnode2loc,basisnode3loc,testlocation,basisarea,analytical1R,analyticalrhoR,analyticalgrad1R,analyticalgradR);
                                    integralA=complex0;
                                    integralB[0]=complex0;integralB[1]=complex0;integralB[2]=complex0;
                                    integralC[0]=complex0;integralC[1]=complex0;integralC[2]=complex0;
                                    for(int indb=0;indb<num_quadpnts;indb++){
                                        basisweight=basisweights[indb];
                                        vec_R[0]=testlocation[0]-basislocations[indb][0];
                                        vec_R[1]=testlocation[1]-basislocations[indb][1];
                                        vec_R[2]=testlocation[2]-basislocations[indb][2];
                                        distance=sqrt(vec_R[0]*vec_R[0]+vec_R[1]*vec_R[1]+vec_R[2]*vec_R[2]);
                                        if(distance>1e-15){
                                            extractedgreensAB=(exp(i*val_kO*distance)-1.0)/(4*pi*distance);
                                            double kRto2=(val_kO*distance)*(val_kO*distance);
                                            double Rto3=distance*distance*distance;
                                            extractedgreensC=((1.0-i*val_kO*distance)*exp(i*val_kO*distance)-1.0-0.5*kRto2)/(4.0*pi*Rto3);
                                        }
                                        else{
                                             extractedgreensAB=i*val_kO/(4.0*pi);
                                             extractedgreensC=0.0;
                                        }
                                        integralA=integralA+basisweight*extractedgreensAB;
                                        double val_dot1;
                                        double vec_val1[3], vec_proj1[3];
                                        vec_val1[0]=-vec_R[0];
                                        vec_val1[1]=-vec_R[1];
                                        vec_val1[2]=-vec_R[2];
                                        val_dot1=dotProduct3d(basistrinormal,vec_val1);
                                        vec_proj1[0]=vec_val1[0]-val_dot1*basistrinormal[0];
                                        vec_proj1[1]=vec_val1[1]-val_dot1*basistrinormal[1];
                                        vec_proj1[2]=vec_val1[2]-val_dot1*basistrinormal[2];
                                        integralB[0]=integralB[0]+vec_proj1[0]*basisweight*extractedgreensAB;
                                        integralB[1]=integralB[1]+vec_proj1[1]*basisweight*extractedgreensAB;
                                        integralB[2]=integralB[2]+vec_proj1[2]*basisweight*extractedgreensAB;
                                        integralC[0]=integralC[0]+basisweight*extractedgreensC*vec_R[0];
                                        integralC[1]=integralC[1]+basisweight*extractedgreensC*vec_R[1];
                                        integralC[2]=integralC[2]+basisweight*extractedgreensC*vec_R[2];
                                    }
                                    integralA+=analytical1R[0];
                                    integralB[0]+=analyticalrhoR[0];
                                    integralB[1]+=analyticalrhoR[1];
                                    integralB[2]+=analyticalrhoR[2];
                                    integralC[0]+=analyticalgrad1R[0]-(val_kO*val_kO)*analyticalgradR[0]/2.0;
                                    integralC[1]+=analyticalgrad1R[1]-(val_kO*val_kO)*analyticalgradR[1]/2.0;
                                    integralC[2]+=analyticalgrad1R[2]-(val_kO*val_kO)*analyticalgradR[2]/2.0;
                                    for(int ind1=0;ind1<3;ind1++){
                                        testfunction=tri2rwgindex[testtri*3+ind1];
                                        if(testfunction!=-1){
                                            if(rwg2box[testfunction]==box_test){
                                                testnode=tri2rwgnode[testtri*3+ind1];
                                                testsign=tri2rwgsign[testtri*3+ind1];
                                                testlength=rwg_lengths[testfunction];
                                                test_node_pos[0]=node_pos[testnode*3+0];
                                                test_node_pos[1]=node_pos[testnode*3+1];
                                                test_node_pos[2]=node_pos[testnode*3+2];
                                                vectorvalue2[0]=testlocation[0]-test_node_pos[0];
                                                vectorvalue2[1]=testlocation[1]-test_node_pos[1];
                                                vectorvalue2[2]=testlocation[2]-test_node_pos[2];
                                                val_dot1=dotProduct3d(basistrinormal,vectorvalue2);
                                                projectionvalue2[0]=vectorvalue2[0]-val_dot1*basistrinormal[0];
                                                projectionvalue2[1]=vectorvalue2[1]-val_dot1*basistrinormal[1];
                                                projectionvalue2[2]=vectorvalue2[2]-val_dot1*basistrinormal[2];
                                                for(int ind2=0;ind2<3;ind2++){
                                                    basisfunction=tri2rwgindex[basistri*3+ind2];
                                                    if(basisfunction!=-1){
                                                        if(rwg2box[basisfunction]==box_basis){
                                                            basisnode=tri2rwgnode[basistri*3+ind2];
                                                            basissign=tri2rwgsign[basistri*3+ind2];
                                                            basislength=rwg_lengths[basisfunction];
                                                            basis_node_pos[0]=node_pos[basisnode*3+0];
                                                            basis_node_pos[1]=node_pos[basisnode*3+1];
                                                            basis_node_pos[2]=node_pos[basisnode*3+2];

                                                            vectorvalue3[0]=testlocation[0]-basis_node_pos[0];
                                                            vectorvalue3[1]=testlocation[1]-basis_node_pos[1];
                                                            vectorvalue3[2]=testlocation[2]-basis_node_pos[2];

                                                            // CALCULATE T-tan OPERATOR
                                                            val_dot2=dotProduct3d(basistrinormal,vectorvalue3);
                                                            projectionvalue3[0]=vectorvalue3[0]-val_dot2*basistrinormal[0];
                                                            projectionvalue3[1]=vectorvalue3[1]-val_dot2*basistrinormal[1];
                                                            projectionvalue3[2]=vectorvalue3[2]-val_dot2*basistrinormal[2];
                                                            double val_efie1=testsign*basissign*testlength*basislength;
                                                            complex<double> val_efie2=(i*val_kO);
                                                            complex<double> val_efie3=(i*val_kO)/4.0;
                                                            complex<double> val_dot_efie1=projectionvalue2[0]*integralB[0]+projectionvalue2[1]*integralB[1]+projectionvalue2[2]*integralB[2];
                                                            double val_dot_efie2=projectionvalue2[0]*projectionvalue3[0]+projectionvalue2[1]*projectionvalue3[1]+projectionvalue2[2]*projectionvalue3[2];
                                                            complex<double> matrix_Znear_1=val_efie1*integralA/val_efie2;
                                                            complex<double> matrix_Znear_2=val_efie1*val_dot_efie1*val_efie3;
                                                            complex<double> matrix_Znear_3=val_efie1*val_dot_efie2*integralA*val_efie3;
                                                            complex<double> val_T_oper_tan = matrix_Znear_1+matrix_Znear_2+matrix_Znear_3;


                                                            double val_mfie1=testsign*basissign*testlength*basislength;
                                                            complex<double> val_identity_tan=complex0;
                                                            complex<double> val_identity_nor=complex0;
                                                            complex<double> val_Kpv_oper_tan=complex0;
                                                            if(testtri==basistri){
                                                                // CALCULATE I-t OPERATOR
                                                                double val_dot_iden_1=vectorvalue2[0]*vectorvalue3[0]+vectorvalue2[1]*vectorvalue3[1]+vectorvalue2[2]*vectorvalue3[2];
                                                                val_identity_tan=val_mfie1*val_dot_iden_1/(basisarea*4.0);

                                                                // CALCULATE I-n OPERATOR
                                                                double vec_cross_1[3];
                                                                vec_cross_1[0]=testtrinormal[1]*vectorvalue3[2]-testtrinormal[2]*vectorvalue3[1];
                                                                vec_cross_1[1]=testtrinormal[2]*vectorvalue3[0]-testtrinormal[0]*vectorvalue3[2];
                                                                vec_cross_1[2]=testtrinormal[0]*vectorvalue3[1]-testtrinormal[1]*vectorvalue3[0];
                                                                double val_dot_iden_2=vectorvalue2[0]*vec_cross_1[0]+vectorvalue2[1]*vec_cross_1[1]+vectorvalue2[2]*vec_cross_1[2];
                                                                val_identity_nor=val_mfie1*val_dot_iden_2/(basisarea*4.0);
                                                            }
                                                            else{
                                                                // CALCULATE Kpv-tan OPERATOR
                                                                complex<double> vec_cross_mfie2[3];
                                                                vec_cross_mfie2[0]=vectorvalue3[1]*integralC[2]-vectorvalue3[2]*integralC[1];
                                                                vec_cross_mfie2[1]=vectorvalue3[2]*integralC[0]-vectorvalue3[0]*integralC[2];
                                                                vec_cross_mfie2[2]=vectorvalue3[0]*integralC[1]-vectorvalue3[1]*integralC[0];
                                                                complex<double> val_dot_mfie2=vectorvalue2[0]*vec_cross_mfie2[0]+vectorvalue2[1]*vec_cross_mfie2[1]+vectorvalue2[2]*vec_cross_mfie2[2];
                                                                val_Kpv_oper_tan=val_mfie1*val_dot_mfie2/4.0;
                                                            }

                                                            //printf("val_imp[%d]: %1.15f %1.15f \n",indbasistri, val_imp.real(), val_imp.imag());

                                                            Znear_index=nearindex[(box_test-box_start)*27+ind_nearbox];
                                                            Znear_index+=(rwg_reverseind[testfunction]-children[box_test*2+0])*num_children[box_basis];
                                                            Znear_index+=rwg_reverseind[basisfunction]-children[box_basis*2+0];

                                                            matrix_Znear11[Znear_index]+=val_T_oper_tan-sqrt(val_imp)*val_identity_tan;
                                                            matrix_Znear12[Znear_index]+=-val_Kpv_oper_tan+(sqrt(1.0/val_imp)-0.5)*val_identity_nor;
                                                            matrix_Znear21[Znear_index]+=val_Kpv_oper_tan+(0.5-sqrt(val_imp))*val_identity_nor;
                                                            matrix_Znear22[Znear_index]+=val_T_oper_tan-sqrt(1.0/val_imp)*val_identity_tan;

/*
                                                            printf("matrix_Znear11[%d]: %1.15f %1.15f \n",Znear_index, matrix_Znear11[Znear_index].real(), matrix_Znear11[Znear_index].imag());
                                                            printf("matrix_Znear12[%d]: %1.15f %1.15f \n",Znear_index, matrix_Znear12[Znear_index].real(), matrix_Znear12[Znear_index].imag());
                                                            printf("matrix_Znear21[%d]: %1.15f %1.15f \n",Znear_index, matrix_Znear21[Znear_index].real(), matrix_Znear21[Znear_index].imag());
                                                            printf("matrix_Znear22[%d]: %1.15f %1.15f \n",Znear_index, matrix_Znear22[Znear_index].real(), matrix_Znear22[Znear_index].imag());
                                                            printf("\n");
                                                            printf("val_T_oper_tan[%d]: %1.15f %1.15f \n",Znear_index, val_T_oper_tan.real(), val_T_oper_tan.imag());
                                                            printf("val_identity_tan[%d]: %1.15f %1.15f \n",Znear_index, val_identity_tan.real(), val_identity_tan.imag());
                                                            printf("sqrt(val_imp)[%d]: %1.15f %1.15f \n",Znear_index, real(sqrt(val_imp)), imag(sqrt(val_imp)));
                                                            printf("\n");
*/

                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            timer_total+=omp_get_wtime()-timer_in;
        }
        double timer_end = omp_get_wtime();
//        fprintf(file_timer,"Znear %d %d %1.5f %1.5f\n", num_mlfmalevels-1, omp_get_thread_num(), timer_end-timer_start, timer_total);
    }
    
    t2 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t2 - t1);
    printf("Znear Calculation Takes %f s \n", time_span.count());

    std::cout.flush();

}







