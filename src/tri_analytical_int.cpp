#include "tri_analytical_int.h"

extern double pi;

void tri_analytical_int(double P1[],double P2[],double P3[],double Pobs[],double triarea,double intval1[],double intval2[],double intval3[],double intval4[]){

//    double P1[3], P2[3], P3[3], Pobs[3];
    double vec_n[3]={0.0, 0.0, 0.0};
    double vec_n_hat[3]={0.0, 0.0, 0.0};
    double temp_vec1[3],temp_vec2[3],temp_vec3[3];
    double temp_val1,temp_val2,temp_val3;
    double vec_l1[3],vec_l2[3],vec_l3[3];
    double vec_l1_hat[3],vec_l2_hat[3],vec_l3_hat[3];
    double vec_u1[3],vec_u2[3],vec_u3[3];
    double vec_u1_hat[3],vec_u2_hat[3],vec_u3_hat[3];
    double vec_r1_Pobs[3],vec_r2_Pobs[3],vec_r3_Pobs[3];
    double z, z_sign;
    double P0_1,P0_2,P0_3;
    double P0_hat_1_dot_u1_hat,P0_hat_2_dot_u2_hat,P0_hat_3_dot_u3_hat;
    double R0_1,R0_2,R0_3;
    double lm_1,lp_1,lm_2,lp_2,lm_3,lp_3;
    double Rm_1,Rp_1,Rm_2,Rp_2,Rm_3,Rp_3;
    double beta_1,beta_2,beta_3;
    double f2_1,f2_2,f2_3;
    double f3_1,f3_2,f3_3;

    intval1[0]=0.0;
    intval2[0]=0.0;
    intval2[1]=0.0;
    intval2[2]=0.0;
    intval3[0]=0.0;
    intval3[1]=0.0;
    intval3[2]=0.0;
    intval4[0]=0.0;
    intval4[1]=0.0;
    intval4[2]=0.0;


//    double intval1=0;
//    double intval2[3]={0.0, 0.0, 0.0};
//    double intval3[3]={0.0, 0.0, 0.0};
//    double intval4[3]={0.0, 0.0, 0.0};
//
//    P1[0]=1.0;
//    P1[1]=0.0;
//    P1[2]=0.0;
//    P2[0]=0.0;
//    P2[1]=2.0;
//    P2[2]=0.0;
//    P3[0]=0.0;
//    P3[1]=0.0;
//    P3[2]=3.0;
//    Pobs[0]=0.5;
//    Pobs[1]=1.0;
//    Pobs[2]=2.0;

    subtract3d(P2,P1,temp_vec1);
    subtract3d(P3,P1,temp_vec2);
    crossProduct3d(temp_vec1, temp_vec2, vec_n);
    normalize3d(vec_n,vec_n_hat);

    subtract3d(P2,P1,vec_l1);
    subtract3d(P3,P2,vec_l2);
    subtract3d(P1,P3,vec_l3);
    normalize3d(vec_l1,vec_l1_hat);
    normalize3d(vec_l2,vec_l2_hat);
    normalize3d(vec_l3,vec_l3_hat);

    crossProduct3d(vec_l1_hat,vec_n_hat,vec_u1);
    crossProduct3d(vec_l2_hat,vec_n_hat,vec_u2);
    crossProduct3d(vec_l3_hat,vec_n_hat,vec_u3);
    normalize3d(vec_u1,vec_u1_hat);
    normalize3d(vec_u2,vec_u2_hat);
    normalize3d(vec_u3,vec_u3_hat);

    subtract3d(P1,Pobs,vec_r1_Pobs);
    subtract3d(P2,Pobs,vec_r2_Pobs);
    subtract3d(P3,Pobs,vec_r3_Pobs);

    subtract3d(Pobs,P1,temp_vec1);
    temp_val1=dotProduct3d(temp_vec1,vec_n_hat);
    z=fabs(temp_val1);
    z_sign=copysign(1.0,temp_val1);

    temp_val1=dotProduct3d(vec_u1_hat,vec_r1_Pobs);
    temp_val2=dotProduct3d(vec_u2_hat,vec_r2_Pobs);
    temp_val3=dotProduct3d(vec_u3_hat,vec_r3_Pobs);
    P0_1=fabs(temp_val1);
    P0_2=fabs(temp_val2);
    P0_3=fabs(temp_val3);
    P0_hat_1_dot_u1_hat=temp_val1/P0_1;
    P0_hat_2_dot_u2_hat=temp_val2/P0_2;
    P0_hat_3_dot_u3_hat=temp_val3/P0_3;

    R0_1=sqrt(z*z + P0_1*P0_1);
    R0_2=sqrt(z*z + P0_2*P0_2);
    R0_3=sqrt(z*z + P0_3*P0_3);

    lm_1=dotProduct3d(vec_r1_Pobs,vec_l1_hat);
    lp_1=dotProduct3d(vec_r2_Pobs,vec_l1_hat);
    lm_2=dotProduct3d(vec_r2_Pobs,vec_l2_hat);
    lp_2=dotProduct3d(vec_r3_Pobs,vec_l2_hat);
    lm_3=dotProduct3d(vec_r3_Pobs,vec_l3_hat);
    lp_3=dotProduct3d(vec_r1_Pobs,vec_l3_hat);

    Rm_1 = norm3d(vec_r1_Pobs);
    Rp_1 = norm3d(vec_r2_Pobs);
    Rm_2 = norm3d(vec_r2_Pobs);
    Rp_2 = norm3d(vec_r3_Pobs);
    Rm_3 = norm3d(vec_r3_Pobs);
    Rp_3 = norm3d(vec_r1_Pobs);

    beta_1 = atan((P0_1*lp_1)/(R0_1*R0_1 + z*Rp_1)) - atan((P0_1*lm_1)/(R0_1*R0_1 + z*Rm_1));
    beta_2 = atan((P0_2*lp_2)/(R0_2*R0_2 + z*Rp_2)) - atan((P0_2*lm_2)/(R0_2*R0_2 + z*Rm_2));
    beta_3 = atan((P0_3*lp_3)/(R0_3*R0_3 + z*Rp_3)) - atan((P0_3*lm_3)/(R0_3*R0_3 + z*Rm_3));

    // Testing the limit cases
    if ((fabs(Rp_1+lp_1)<1e-15)&&(fabs(Rm_1+lm_1)<1e-15))
        f2_1=log(lm_1/lp_1);
    else
        f2_1=log((Rp_1+lp_1)/(Rm_1+lm_1));
    if ((fabs(Rp_2+lp_2)<1e-15)&&(fabs(Rm_2+lm_2)<1e-15))
        f2_2=log(lm_2/lp_2);
    else
        f2_2=log((Rp_2+lp_2)/(Rm_2+lm_2));
    if ((fabs(Rp_3+lp_3)<1e-15)&&(fabs(Rm_3+lm_3)<1e-15))
        f2_3=log(lm_3/lp_3);
    else
        f2_3=log((Rp_3+lp_3)/(Rm_3+lm_3));

    f3_1=R0_1*R0_1*f2_1+lp_1*Rp_1-lm_1*Rm_1;
    f3_2=R0_2*R0_2*f2_2+lp_2*Rp_2-lm_2*Rm_2;
    f3_3=R0_3*R0_3*f2_3+lp_3*Rp_3-lm_3*Rm_3;


    // integral (1/R)
    if (P0_1>0){intval1[0]=intval1[0]+P0_hat_1_dot_u1_hat*(P0_1*f2_1-z*beta_1);}
    if (P0_2>0){intval1[0]=intval1[0]+P0_hat_2_dot_u2_hat*(P0_2*f2_2-z*beta_2);}
    if (P0_3>0){intval1[0]=intval1[0]+P0_hat_3_dot_u3_hat*(P0_3*f2_3-z*beta_3);}

    // integral (rhoprime/R)
    multiply3d(vec_u1_hat,f3_1/2,temp_vec1);
    add3d(intval2,temp_vec1,intval2);
    multiply3d(vec_u2_hat,f3_2/2,temp_vec2);
    add3d(intval2,temp_vec2,intval2);
    multiply3d(vec_u3_hat,f3_3/2,temp_vec3);
    add3d(intval2,temp_vec3,intval2);

    // integral gradprime(1/R)
    if (P0_1>0){
        intval3[0]+=(vec_u1_hat[0]*f2_1)+(vec_n_hat[0]*z_sign*beta_1*P0_hat_1_dot_u1_hat);
        intval3[1]+=(vec_u1_hat[1]*f2_1)+(vec_n_hat[1]*z_sign*beta_1*P0_hat_1_dot_u1_hat);
        intval3[2]+=(vec_u1_hat[2]*f2_1)+(vec_n_hat[2]*z_sign*beta_1*P0_hat_1_dot_u1_hat);
    }
    if (P0_2>0){
        intval3[0]+=(vec_u2_hat[0]*f2_2)+(vec_n_hat[0]*z_sign*beta_2*P0_hat_2_dot_u2_hat);
        intval3[1]+=(vec_u2_hat[1]*f2_2)+(vec_n_hat[1]*z_sign*beta_2*P0_hat_2_dot_u2_hat);
        intval3[2]+=(vec_u2_hat[2]*f2_2)+(vec_n_hat[2]*z_sign*beta_2*P0_hat_2_dot_u2_hat);
    }
    if (P0_3>0){
        intval3[0]+=(vec_u3_hat[0]*f2_3)+(vec_n_hat[0]*z_sign*beta_3*P0_hat_3_dot_u3_hat);
        intval3[1]+=(vec_u3_hat[1]*f2_3)+(vec_n_hat[1]*z_sign*beta_3*P0_hat_3_dot_u3_hat);
        intval3[2]+=(vec_u3_hat[2]*f2_3)+(vec_n_hat[2]*z_sign*beta_3*P0_hat_3_dot_u3_hat);
    }

    // integral gradprime(R)
    if (P0_1>0){
        intval4[0]+=(vec_u1_hat[0]*(f3_1)/2)+(vec_n_hat[0]*z_sign*z*z*P0_hat_1_dot_u1_hat*beta_1)-(vec_n_hat[0]*z_sign*z*P0_hat_1_dot_u1_hat*P0_1*f2_1);
        intval4[1]+=(vec_u1_hat[1]*(f3_1)/2)+(vec_n_hat[1]*z_sign*z*z*P0_hat_1_dot_u1_hat*beta_1)-(vec_n_hat[1]*z_sign*z*P0_hat_1_dot_u1_hat*P0_1*f2_1);
        intval4[2]+=(vec_u1_hat[2]*(f3_1)/2)+(vec_n_hat[2]*z_sign*z*z*P0_hat_1_dot_u1_hat*beta_1)-(vec_n_hat[2]*z_sign*z*P0_hat_1_dot_u1_hat*P0_1*f2_1);
    }
    if (P0_2>0){
        intval4[0]+=(vec_u2_hat[0]*(f3_2)/2)+(vec_n_hat[0]*z_sign*z*z*P0_hat_2_dot_u2_hat*beta_2)-(vec_n_hat[0]*z_sign*z*P0_hat_2_dot_u2_hat*P0_2*f2_2);
        intval4[1]+=(vec_u2_hat[1]*(f3_2)/2)+(vec_n_hat[1]*z_sign*z*z*P0_hat_2_dot_u2_hat*beta_2)-(vec_n_hat[1]*z_sign*z*P0_hat_2_dot_u2_hat*P0_2*f2_2);
        intval4[2]+=(vec_u2_hat[2]*(f3_2)/2)+(vec_n_hat[2]*z_sign*z*z*P0_hat_2_dot_u2_hat*beta_2)-(vec_n_hat[2]*z_sign*z*P0_hat_2_dot_u2_hat*P0_2*f2_2);
    }
    if (P0_3>0){
        intval4[0]+=(vec_u3_hat[0]*(f3_3)/2)+(vec_n_hat[0]*z_sign*z*z*P0_hat_3_dot_u3_hat*beta_3)-(vec_n_hat[0]*z_sign*z*P0_hat_3_dot_u3_hat*P0_3*f2_3);
        intval4[1]+=(vec_u3_hat[1]*(f3_3)/2)+(vec_n_hat[1]*z_sign*z*z*P0_hat_3_dot_u3_hat*beta_3)-(vec_n_hat[1]*z_sign*z*P0_hat_3_dot_u3_hat*P0_3*f2_3);
        intval4[2]+=(vec_u3_hat[2]*(f3_3)/2)+(vec_n_hat[2]*z_sign*z*z*P0_hat_3_dot_u3_hat*beta_3)-(vec_n_hat[2]*z_sign*z*P0_hat_3_dot_u3_hat*P0_3*f2_3);
    }


    intval1[0]/=(4*pi*triarea);
    intval2[0]/=(4*pi*triarea);
    intval2[1]/=(4*pi*triarea);
    intval2[2]/=(4*pi*triarea);
    intval3[0]/=(4*pi*triarea);
    intval3[1]/=(4*pi*triarea);
    intval3[2]/=(4*pi*triarea);
    intval4[0]/=(4*pi*triarea);
    intval4[1]/=(4*pi*triarea);
    intval4[2]/=(4*pi*triarea);

}
