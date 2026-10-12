#include "calculate_translationoperator.h"

void calculate_translationoperator(){
    //
    //
    // calculate_translationoperator
    // Calculate the translation operator

    int ind_trans_op=0;
    double translationvector[3]{0.0, 0.0, 0.0};
    double translationvectorsize=0.0;
    double hankelarg=0;
    double legendrearg=0;
    complex<double> val_i=complex<double>(0.0, 1.0);
    double *thetalocations;
    double *thetaweights;
    double *philocations;
    double *phiweights;
    double *sintheta;
    double *costheta;
    double *sinphi;
    double *cosphi;
    double khat[3]{0.0, 0.0, 0.0};
    complex<double> *hankelvalue;
    double *legendrevalue;

    //
    // Translation operator for MLFMA (mvm_Zfardotx)
    //
    trans_op=new complex<double>[size_trans];
    for(int ind=0; ind<size_trans; ind++){
        trans_op[ind]=complex<double>(0.0, 0.0);
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    for(int indlev=2; indlev<num_mlfmalevels; indlev++){
        int numberofterms=num_truncation[indlev];
        int numberofthetas=numberofterms+1;
        int numberofphis=numberofthetas*2;

        thetalocations=new double[numberofthetas];
        thetaweights=new double[numberofthetas];
        philocations=new double[numberofphis];
        phiweights=new double[numberofphis];
        sintheta=new double[numberofthetas];
        costheta=new double[numberofthetas];
        sinphi=new double[numberofphis];
        cosphi=new double[numberofphis];
        theta_locations_weights(numberofthetas, thetaweights, thetalocations);
        phi_locations_weights(numberofphis, phiweights, philocations);
        calculate_sin_cos(numberofthetas,numberofphis, thetalocations, philocations, sintheta, costheta, sinphi, cosphi);

        omp_set_num_threads(num_threads);
        #pragma omp parallel if(parallel_omp_on)
        {
            double timer_start = omp_get_wtime();
            double timer_total = 0.0;
            #pragma omp for private(ind_trans_op,translationvector,translationvectorsize,hankelarg,\
                                             hankelvalue,khat,legendrearg,legendrevalue) schedule(dynamic, 1)
            for(int ind=0;ind<(7*7*7);ind++){
                double timer_in = omp_get_wtime();
                hankelvalue=new complex<double>[numberofterms+1];
                legendrevalue=new double[numberofterms+1];
                int indz=ind%7;
                int indy=(ind/7)%7;
                int indx=(ind/49)%7;
                int indexval=indlev*343+ind;
                if(trans_index[indexval]!=-1){
                    ind_trans_op=trans_index[indexval];
                    translationvector[0]=(indx-3)*size_boxes[indlev];
                    translationvector[1]=(indy-3)*size_boxes[indlev];
                    translationvector[2]=(indz-3)*size_boxes[indlev];
                    translationvectorsize=sqrt(translationvector[0]*translationvector[0] \
                                              +translationvector[1]*translationvector[1] \
                                              +translationvector[2]*translationvector[2]);
                    hankelarg=val_kO*translationvectorsize;
                    for(int indterm=0; indterm<(numberofterms+1); indterm++){
                        hankelvalue[indterm].real(sph_bessel(indterm, hankelarg));
                        hankelvalue[indterm].imag(sph_neumann(indterm, hankelarg));
                    }
                    for(int indtheta=0; indtheta<numberofthetas; indtheta++){
                        for(int indphi=0; indphi<numberofphis; indphi++){
                            khat[0]=sintheta[indtheta]*cosphi[indphi];
                            khat[1]=sintheta[indtheta]*sinphi[indphi];
                            khat[2]=costheta[indtheta];
                            legendrearg=(translationvector[0]*khat[0]+translationvector[1]*khat[1]+translationvector[2]*khat[2])/translationvectorsize;
                            legendrevalue[0]=1.0;
                            legendrevalue[1]=legendrearg;
                            for(int indterm=2; indterm<(numberofterms+1); indterm++){
                                legendrevalue[indterm]=((2*indterm-1)*legendrearg*legendrevalue[indterm-1] - (indterm-1)*legendrevalue[indterm-2])/indterm;
    //                                    legendrevalue[indterm]=legendrel(indterm,legendrearg);
                            }
                            for(int indterm=0; indterm<numberofterms+1; indterm++){
                                trans_op[ind_trans_op]+=pow(val_i,indterm)*double(2*indterm+1)*hankelvalue[indterm]*legendrevalue[indterm];
                            }
                            trans_op[ind_trans_op]*=(thetaweights[indtheta]*phiweights[indphi]);
                            ind_trans_op++;
                        }
                    }
                }
                delete[] hankelvalue;
                delete[] legendrevalue;
                timer_total+=omp_get_wtime()-timer_in;
            }
            double timer_end = omp_get_wtime();
//            fprintf(file_timer,"Trans%d %d %1.5f %1.5f \n", indlev, omp_get_thread_num(), timer_end-timer_start, timer_total);
        }
        delete[] thetalocations;
        delete[] thetaweights;
        delete[] philocations;
        delete[] phiweights;
        delete[] sintheta;
        delete[] costheta;
        delete[] sinphi;
        delete[] cosphi;

    }
    printf("Translation Operator Size: %d  [%f GB]\n", size_trans, size_trans*16.0/1024.0/1024/1024);


    //
    // Translation operator for AMLFMA Preconditioner (mvm_Zfardotx_in)
    //
    if(val_pc==2 && val_approx>0){
        trans_op_pc=new complex<double>[size_trans_pc];
        for(int ind=0; ind<size_trans_pc; ind++){
            trans_op_pc[ind]=complex<double>(0.0, 0.0);
        }

        for(int indlev=2; indlev<num_mlfmalevels; indlev++){
            int numberofterms=num_truncation_pc[indlev];
            int numberofthetas=numberofterms+1;
            int numberofphis=numberofthetas*2;

            thetalocations=new double[numberofthetas];
            thetaweights=new double[numberofthetas];
            philocations=new double[numberofphis];
            phiweights=new double[numberofphis];
            sintheta=new double[numberofthetas];
            costheta=new double[numberofthetas];
            sinphi=new double[numberofphis];
            cosphi=new double[numberofphis];
            theta_locations_weights(numberofthetas, thetaweights, thetalocations);
            phi_locations_weights(numberofphis, phiweights, philocations);
            calculate_sin_cos(numberofthetas,numberofphis, thetalocations, philocations, sintheta, costheta, sinphi, cosphi);

            hankelvalue=new complex<double>[numberofterms+1];
            legendrevalue=new double[numberofterms+1];
            for(int ind=0; ind<(numberofterms+1); ind++){
                legendrevalue[ind]=0.0;
                hankelvalue[ind]=complex<double>(0.0, 0.0);
            }

            omp_set_num_threads(num_threads);
            #pragma omp parallel if(parallel_omp_on)
            {
                double timer_start = omp_get_wtime();
                double timer_total = 0.0;
                #pragma omp for private(ind_trans_op,translationvector,translationvectorsize,hankelarg,\
                                                 hankelvalue,khat,legendrearg,legendrevalue) schedule(dynamic, 1)

                for(int ind=0;ind<(7*7*7);ind++){
                    double timer_in = omp_get_wtime();
                    hankelvalue=new complex<double>[numberofterms+1];
                    legendrevalue=new double[numberofterms+1];
                    int indz=ind%7;
                    int indy=(ind/7)%7;
                    int indx=(ind/49)%7;
                    int indexval=indlev*343+ind;
                    if(trans_index_pc[indexval]!=-1){
                        ind_trans_op=trans_index_pc[indexval];
                        translationvector[0]=(indx-3)*size_boxes[indlev];
                        translationvector[1]=(indy-3)*size_boxes[indlev];
                        translationvector[2]=(indz-3)*size_boxes[indlev];
                        translationvectorsize=sqrt(translationvector[0]*translationvector[0] \
                                                  +translationvector[1]*translationvector[1] \
                                                  +translationvector[2]*translationvector[2]);
                        hankelarg=val_kO*translationvectorsize;
                        for(int indterm=0; indterm<(numberofterms+1); indterm++){
                            hankelvalue[indterm].real(sph_bessel(indterm, hankelarg));
                            hankelvalue[indterm].imag(sph_neumann(indterm, hankelarg));
                        }
                        for(int indtheta=0; indtheta<numberofthetas; indtheta++){
                            for(int indphi=0; indphi<numberofphis; indphi++){
                                khat[0]=sintheta[indtheta]*cosphi[indphi];
                                khat[1]=sintheta[indtheta]*sinphi[indphi];
                                khat[2]=costheta[indtheta];
                                legendrearg=(translationvector[0]*khat[0]+translationvector[1]*khat[1]+translationvector[2]*khat[2])/translationvectorsize;
                                legendrevalue[0]=1.0;
                                legendrevalue[1]=legendrearg;
                                for(int indterm=2; indterm<(numberofterms+1); indterm++){
                                    legendrevalue[indterm]=((2*indterm-1)*legendrearg*legendrevalue[indterm-1] - (indterm-1)*legendrevalue[indterm-2])/indterm;
    //                                    legendrevalue[indterm]=legendrel(indterm,legendrearg);
                                }
                                for(int indterm=0; indterm<numberofterms+1; indterm++){
                                    trans_op_pc[ind_trans_op]+=pow(val_i,indterm)*double(2*indterm+1)*hankelvalue[indterm]*legendrevalue[indterm];
                                }
                                trans_op_pc[ind_trans_op]*=(thetaweights[indtheta]*phiweights[indphi]);
                                ind_trans_op++;
                            }
                        }
                    }
                    delete[] hankelvalue;
                    delete[] legendrevalue;
                    timer_total+=omp_get_wtime()-timer_in;
                }
                double timer_end = omp_get_wtime();
//                fprintf(file_timer,"TransPC%d %d %1.5f %1.5f\n", indlev, omp_get_thread_num(), timer_end-timer_start, timer_total);
            }
            delete[] thetalocations;
            delete[] thetaweights;
            delete[] philocations;
            delete[] phiweights;
            delete[] sintheta;
            delete[] costheta;
            delete[] sinphi;
            delete[] cosphi;
        }
//        printf("PC Translation Operator Size: %d  [%f GB]\n", size_trans_pc, size_trans_pc*16.0/1024.0/1024/1024);
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    printf("Translation Operator Calculation Takes %f s \n", time_span.count());
//    fclose(file_timer);


    std::cout.flush();
}

