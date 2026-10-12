#include <complex>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <iostream>

using namespace std;

extern char modelname[];
extern int num_mlfmalevels;
extern int num_digits;
extern int num_interpol;
extern int source_type;
extern int val_itermax;
extern double val_soltol;
extern int num_media;
extern double val_cfiealpha[];
extern complex<double> val_impedance[];
extern double val_freq;
extern int num_freqs;
extern double freq_start, freq_end;
extern complex<double> *source_CSP;
extern complex<double> *source_current;
extern complex<double> *source_DG;
extern double *source_dipole;
extern double source_PW[];
extern int *num_source;
extern int val_pol;
extern int val_pc;
extern double val_approx;
extern int val_itermax_pc;
extern int val_typeff;
extern int *num_ffdir_tp;
extern double dir_ff;
extern double dir_cfc[];
extern int parallel_omp_on;
extern int num_threads;

extern int num_nodes;
extern int num_tris;
extern double *node_pos;
extern int *tri_nodes;
extern int *tri_media;

extern double pi;
extern double val_e;
extern double val_epsilonO;
extern double val_muO;
extern double val_omega;
extern double val_c;
extern double val_lambdaO;
extern double val_kO;
extern double val_etaO;

void read_simulationinput(){
    //
    //
    // read_simulationinput
    // Setting the input parameters for the simulation


    int filestat;
    char char_temp[200];
    FILE *pFile;
    int ind_formulation;
    double val_real, val_imag;

    pFile = fopen ("input.pnu","r");
    if (pFile==NULL){
        printf("Input File Not Found! \n");
        exit(EXIT_FAILURE);
    }

    //Reading until the simulation inputs are listed
    while(!feof(pFile)){
        filestat=fscanf(pFile, "%s\n", char_temp); // flag: "SimulationInputs"
        if(!strcmp(char_temp,"SimulationInputs")){
            break;
        }
    }

    //
    //Simulation parameters are listed. 
    //These parameters are read in line order.
    filestat=fscanf(pFile, "%s %s\n", char_temp, modelname);
    printf("model_name %s \n", modelname);
    filestat=fscanf(pFile, "%s %d\n", char_temp, &ind_formulation);
    printf("Formulations %d \n", ind_formulation);
    filestat=fscanf(pFile, "%s %lf %lf %d\n", char_temp, &freq_start, &freq_end, &num_freqs);
    printf("Frequency %lf %lf %d \n", freq_start, freq_end, num_freqs);
    filestat=fscanf(pFile, "%s %d\n", char_temp, &source_type);
    printf("source_type %d \n", source_type);
    if(source_type==0){
        filestat=fscanf(pFile, "%s %lf %lf %d %lf %lf %d %d\n", char_temp, &source_PW[0], &source_PW[1], &num_source[0], &source_PW[2], &source_PW[3], &num_source[1], &val_pol);
        printf("source_PW %lf %lf %d %lf %lf %d %d \n", source_PW[0], source_PW[1], num_source[0], source_PW[2], source_PW[3], num_source[1], val_pol);
    } else{
        printf("Error: wrong excitation input! \n");
        exit(EXIT_FAILURE);
    }
    filestat=fscanf(pFile, "%s %d\n", char_temp, &val_pc);
    printf("Preconditioners %d \n", val_pc);
    filestat=fscanf(pFile, "%s %d\n", char_temp, &num_mlfmalevels);
    printf("num_mlfmalevels %d \n", num_mlfmalevels);
    filestat=fscanf(pFile, "%s %lf\n", char_temp, &dir_ff);
    printf("dir_ff %lf \n", dir_ff);
    filestat=fscanf(pFile, "%s %d\n", char_temp, &num_digits);
    printf("num_digits %d \n", num_digits);
    filestat=fscanf(pFile, "%s %lf\n", char_temp, &val_approx);
    printf("val_approx %lf \n", val_approx);
    filestat=fscanf(pFile, "%s %d\n", char_temp, &num_interpol);
    printf("num_interpol %d \n", num_interpol);
    filestat=fscanf(pFile, "%s %d\n", char_temp, &val_itermax_pc);
    printf("val_itermax_pc %d \n", val_itermax_pc);
    filestat=fscanf(pFile, "%s %lf\n", char_temp, &val_soltol);
    printf("val_approx %lf \n", val_soltol);
    filestat=fscanf(pFile, "%s %lf %lf\n", char_temp, &dir_cfc[0], &dir_cfc[1]);
    printf("dir_cfc %lf %lf\n", dir_cfc[0], dir_cfc[1]);
    filestat=fscanf(pFile, "%s %d\n", char_temp, &val_itermax);
    printf("val_itermax %d \n", val_itermax);
    filestat=fscanf(pFile, "%s %d %d\n", char_temp, &parallel_omp_on, &num_threads);
    printf("parallel_omp %d %d\n",  parallel_omp_on, num_threads);
    filestat=fscanf(pFile, "%s %s\n", char_temp, char_temp); // flag: "## End-Param"
    filestat=fscanf(pFile, "%s %d\n", char_temp, &num_media);
    printf("ModelParts %d \n", num_media);
    if(ind_formulation==0){
        for (int i=0; i<num_media; i++){
            filestat=fscanf(pFile, "%s %lf\n", char_temp, &val_cfiealpha[i]); 
            printf("val_cfiealpha[%d]: %f \n", i, val_cfiealpha[i]);
        }
    }
    else if(ind_formulation==1){
        for (int i=0; i<num_media; i++){
            filestat=fscanf(pFile, "%s %lf %lf\n", char_temp, &val_real, &val_imag); 
            val_impedance[i].real(val_real);
            val_impedance[i].imag(val_imag);
            printf("val_impedance[%d]: %f %f \n", i, val_impedance[i].real(), val_impedance[i].imag());
        }
    }

    //
    //Mesh model details are read.
    //The points are cells are read, respectively.
    //
    //Reading point coordinates
    while(!feof(pFile)){
        filestat=fscanf(pFile, "%s", char_temp);
        if(!strcmp(char_temp,"POINTS")){
            break;
        }
    }
    filestat=fscanf(pFile, "%d \n", &num_nodes);
    printf("Number of points: %d \n", num_nodes);
    node_pos = new double[num_nodes*3];
    for(int indline=0;indline<num_nodes;indline++){
        filestat=fscanf (pFile, "%lf %lf %lf\n", &node_pos[indline*3+0], &node_pos[indline*3+1], &node_pos[indline*3+2]);
    }
    //
    // Reading cell points
    while(!feof(pFile)){
        filestat=fscanf(pFile, "%s", char_temp);
        if(!strcmp(char_temp,"CELLS")){
            break;
        }
    }
    filestat=fscanf(pFile, "%d \n", &num_tris);
    printf("Number of cells: %d \n", num_tris);
    tri_nodes = new int[num_tris*3];
    tri_media = new int[num_tris];
    for(int indline=0;indline<num_tris;indline++){
        filestat=fscanf(pFile, "%d %d %d %d\n", &tri_media[indline], &tri_nodes[indline*3+0], &tri_nodes[indline*3+1], &tri_nodes[indline*3+2]);
        //tri_media[indline]--;
    }
    //
    fclose (pFile);





    ////////
    ///////
    //////
    /////
    ////
    ///
    //
    val_freq=freq_start;
    //
    ///
    ////
    /////
    //////
    ///////
    ////////


    //Generating the global constants
    val_omega = 2*pi*val_freq;
    val_c = 1/sqrt(val_epsilonO*val_muO);
    val_lambdaO = val_c/val_freq;
    val_kO = 2*pi/val_lambdaO;
    val_etaO = sqrt(val_muO/val_epsilonO);


    printf("Frequency: %f GHz\n", val_freq/1e9);
    printf("Illumination Type: ");
    if(source_type==0){
        printf("Plane Wave \n");
    }else if(source_type==1){
        printf("Delta Gap \n");
    }else if(source_type==2){
        printf("Dipole \n");
    }else if(source_type==3){
        printf("Current Source \n");
    }else if(source_type==4){
        printf("Complex Source Point \n");
    }else{
        printf("Invalid!!! \n");
        exit(EXIT_FAILURE);
    }
    printf("\n");
    printf("PI: %1.20f \n", pi);
    printf("Omega: %f rad Hz \n", val_omega);
    printf("Epsilon: %1.4e F/m \n", val_epsilonO);
    printf("Mu: %1.4e H/m \n", val_muO);
    printf("Speed of Light: %1.4e m/s \n", val_c);
    printf("Wavelength %f m \n", val_lambdaO);
    printf("Wave number: %f rad/m \n", val_kO);
    printf("Wave Impedance: %f ohm \n", val_etaO);
    printf("\n");

    std::cout.flush();



}



