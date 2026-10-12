#include <cmath>
#include <complex>

using namespace std;

void subtract3d(double vec1[], double vec2[], double vec3[]){
    // Subtract given three dimensional vectors vec1 and vec2
    // vec3=vec1-vec2;

    vec3[0]=vec1[0]-vec2[0];
    vec3[1]=vec1[1]-vec2[1];
    vec3[2]=vec1[2]-vec2[2];
}

void add3d(double vec1[], double vec2[], double vec3[]){
    // Add given three dimensional vectors vec1 and vec2
    // vec3=vec1+vec2;

    vec3[0]=vec1[0]+vec2[0];
    vec3[1]=vec1[1]+vec2[1];
    vec3[2]=vec1[2]+vec2[2];
}

void multiply3d(double vec1[], double val1, double vec2[]){
    // Multiply given three dimensional vector vec1 and scalar val1
    // vec2=vec1*val1;

    vec2[0]=vec1[0]*val1;
    vec2[1]=vec1[1]*val1;
    vec2[2]=vec1[2]*val1;
}

void crossProduct3d(double vec1[], double vec2[], double resvec[]){
    // Cross Product of three dimensional vectors, vec1 and vec2.
    // resvec = vec1 x vec2;

    resvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    resvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    resvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

}

void crossProduct3d_cd(complex<double> vec1[], double vec2[], complex<double> resvec[]){
    // Cross Product of three dimensional vectors, vec1 and vec2.
    // resvec = vec1 x vec2;

    resvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    resvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    resvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

}

void crossProduct3d_dc(double vec1[], complex<double> vec2[], complex<double> resvec[]){
    // Cross Product of three dimensional vectors, vec1 and vec2.
    // resvec = vec1 x vec2;

    resvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    resvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    resvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

}

void crossProduct3d_c(complex<double> vec1[], complex<double> vec2[], complex<double> resvec[]){
    // Cross Product of three dimensional vectors, vec1 and vec2.
    // resvec = vec1 x vec2;

    resvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    resvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    resvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

}

double dotProduct3d(double vec1[], double vec2[]){
    // Dot Product of three dimensional vectors, vec1 and vec2
    // resval=vec1.vec2;

    double resval=0;
    resval=vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
    return resval;
}

complex<double> dotProduct3d_dc(double vec1[], complex<double> vec2[]){
    // Dot Product of three dimensional vectors, vec1 and vec2
    // resval=vec1.vec2;

    complex<double> resval=complex<double>{0.0,0.0};
    resval=vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
    return resval;
}

complex<double> dotProduct3d_c(complex<double> vec1[], complex<double> vec2[]){
    // Dot Product of three dimensional vectors, vec1 and vec2
    // resval=vec1.vec2;

    complex<double> resval=complex<double>{0.0,0.0};
    resval=vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
    return resval;
}

void normalize3d(double vec1[], double vec2[]){
    // Normalize a given three dimensional vector vec1
    // vec2=vec1/norm(vec1);

    double val_norm;
    val_norm=sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
    vec2[0]=vec1[0]/val_norm;
    vec2[1]=vec1[1]/val_norm;
    vec2[2]=vec1[2]/val_norm;
}

double norm3d(double vec1[]){
    // Calculate two-norm of a given three dimensional vector vec1
    // norm(vec1);

    double val_norm;
    val_norm=sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
    return val_norm;
}

//double normNd_c(complex<double> vec1[],int valDim){
//    // Calculate two-norm of a given N-dimensional vector vec1
//    // norm(vec1,2);
//
//    double val_norm=0.0;
//    for(int ind=0;ind<valDim;ind++){
//        val_norm+=vec1[ind].real()*vec1[ind].real()+vec1[ind].imag()*vec1[ind].imag();
//    }
//    val_norm=sqrt(val_norm);
//    return val_norm;
//}
