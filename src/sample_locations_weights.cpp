#include <cmath>

extern double pi;

void phi_locations_weights(int num_sample_phi, double *weights_phi, double *locs_phi){
    for(int indphi=0; indphi<num_sample_phi; indphi++){
        weights_phi[indphi]=2*pi/num_sample_phi;
        locs_phi[indphi]=2*pi/num_sample_phi*indphi;
    }
}

void theta_locations_weights(int num_sample_theta, double *weights_theta, double *locs_theta){
    // theta_locations_weights
    // Gauss-Legendre Quadrature Points and Weights
    int s=num_sample_theta;
    int m=(s+1)/2;
    double P1=0.0;
    double P2=0.0;
    double P3=0.0;
    double Pder=0.0;
    double z=0.0;
    double z1=1.0;

    for(int ind_sample=0; ind_sample<m; ind_sample++){
        z=cos(pi*(ind_sample-0.25+1)/(s+0.5)); // First guess of the ind_sample'th root of the s-order Legendre polynomial
        z1=1.0;
        while(std::abs(z1)>1e-15){ // "std::abs" returns float output, "abs" returns int output
            P1=1.0; // Legendre polynomial value
            P2=0.0;
            for(int ind_leg=1; ind_leg<s+1; ind_leg++){
                P3=P2;
                P2=P1;
                P1=((2*ind_leg-1)*z*P2 - (ind_leg-1)*P3)/ind_leg;
            }
            Pder = s*(z*P1-P2)/(z*z-1); // Derivative of Legendre polynomial
            z1=P1/Pder;
            z=z-z1; // Newton-Raphson iteration
        }
        locs_theta[ind_sample]=acos(z);
        locs_theta[s-1-ind_sample]=acos(-z);
        weights_theta[ind_sample]=2/((1- (z*z) )* (Pder*Pder) );
        weights_theta[s-1-ind_sample]=weights_theta[ind_sample];

    }
}

void theta_locations_weights_lin(int num_sample_theta, double *weights_theta, double *locs_theta){
    // theta_locations_weights_lin
    // Linearly Sampled Points and Weights
    for(int indtheta=0; indtheta<num_sample_theta; indtheta++){
        weights_theta[indtheta]=pi/num_sample_theta;
        locs_theta[indtheta]=pi/(num_sample_theta-1)*indtheta;
    }
}

void calculate_sin_cos(int num_thetas, int num_phis, double *thetas, double *phis, double *sintheta, double *costheta, double *sinphi, double *cosphi){

    for(int indtheta=0; indtheta<num_thetas; indtheta++){
        sintheta[indtheta]=sin(thetas[indtheta]);
        costheta[indtheta]=cos(thetas[indtheta]);
    }
    for(int indphi=0; indphi<num_phis; indphi++){
        sinphi[indphi]=sin(phis[indphi]);
        cosphi[indphi]=cos(phis[indphi]);
    }

}



