#include <complex>
#include <cmath>
#include <iostream>
#include <stdio.h>

std::complex<double>** matrix_inverse(std::complex<double>** matrix, int size){
    int i,j,k;
    int n=size;
    std::complex<double> d;


    std::complex<double>** result = new std::complex<double>*[size];
     for(int i = 0; i < size; i++)
     {
         result[i] = new std::complex<double>[size];
         result[i][i]=std::complex<double> (1.0, 0.0);
     }

    /************** partial pivoting **************/
    for(i=(n-1);i>0;i--)
    {
        if(abs(matrix[i-1][0]) < abs(matrix[i][0]))
        for(j=0;j<n;j++)
        {
            d=matrix[i][j];
            matrix[i][j]=matrix[i-1][j];
            matrix[i-1][j]=d;

            // added
            d=result[i][j];
            result[i][j]=result[i-1][j];
            result[i-1][j]=d;
        }
    }


    /********** reducing to diagonal  matrix ***********/

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        if(j!=i)
        {
            d=matrix[j][i]/matrix[i][i];
            for(k=0;k<n;k++)
                matrix[j][k]-=matrix[i][k]*d;

            // added
            for(k=0;k<n;k++)
                result[j][k]-=result[i][k]*d;
        }
    }
    /************** reducing to unit matrix *************/
    for(i=0;i<n;i++)
    {
    d=matrix[i][i];
        for(j=0;j<n;j++)
            matrix[i][j]=matrix[i][j]/d;

        // added
        for(j=0;j<n;j++)
            result[i][j]=result[i][j]/d;
    }

/*
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            printf("Matinv[%d][%d]: %1.16f %1.16f \n", i, j, result[i][j].real(), result[i][j].imag());
*/
    

    return result;
}