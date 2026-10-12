#include <iostream>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <string.h>

 typedef struct
{
	double **M;
	int n, m;
} mmatrix;


void copyAndSetPartOfArray2(double** dst, double** src, int rowSize, int colSize, int dstStartRow, int dstStartCol,  int srcStartRow, int srcStartCol)
{
    // dst: 2d dynamic matrix ; src 2d dynamic matrix
    // Example:
    // F(8:10,1:4) = out(8:10,1:4);
    // copyAndSetPartOfArray2(F, out, 3, 4, 7, 0,  7, 0);
    int i;
    for(i = 0; i<rowSize; i++)
    {
        memcpy(&dst[dstStartRow][dstStartCol], &src[srcStartRow][srcStartCol], sizeof(double)*colSize);
        dstStartRow++;
        srcStartRow++;
    }
}

void freematr(mmatrix *A)
{
	if (A)
	{
		if (A->M != 0)
		{
			for (int i = 0; i < A->n; ++i)
				free(A->M[i]);
			free(A->M);
		}
		free(A);
	}
}

mmatrix *memmatr(int filas, int columnas)
{
	int i;
	mmatrix *A;
	A = (mmatrix*)malloc(sizeof(mmatrix));
	if (!A) return 0;
	A->n = filas;
	A->m = columnas;
	A->M = (double**)malloc(filas*sizeof(double*));
	if (A->M == 0) return 0;
	for (i = 0; i < filas; ++i)
	{
		A->M[i] = (double*)malloc(columnas*sizeof(double));
		if (A->M[i] == 0)
		{
			A->n = i;
			freematr(A);
			return 0;
		}
	}
	return A;
}

void cofacmmatrix(mmatrix *A, mmatrix *D, int p, int q, int tam)
{
	int i, j, k = 0, w = 0;

	for (i = 0; i < tam; ++i){
		for (j = 0; j < tam; ++j){
			if (i != p && j != q)
			{
				D->M[k][w++] = A->M[i][j];
				if (w == (tam-1))
				{
					w = 0;
					k++;
				}
			}
		}
	}
}

double determmatrix(mmatrix *A, int tam)
{
	int i, j, signo = 1;
	double determ=0.0;
	mmatrix *D;

	if(A->n == A->m)
	{
		if (tam == 1) return A->M[0][0];
		if (tam == 2) return (A->M[0][0]* A->M[1][1]- A->M[0][1]* A->M[1][0]);
		else
		{
			D = memmatr(A->n, A->m);
			for (i = 0; i < tam; ++i)
			{
				cofacmmatrix(A, D, 0, i, tam);
				determ += signo*A->M[0][i]*determmatrix(D, (tam - 1));

                //printf("determ: %1.16f \n", determ);
				signo = -signo;
			}
			freematr(D);
		}
	}
    return determ;
}

void adjmmatrix(mmatrix *A, mmatrix *Adj, int tam)
{
	mmatrix *E;
	int i, j, signo2 = 1;

	if (tam == 1) Adj->M[0][0] = 1;

	for (i = 0; i < tam; ++i)
		for (j = 0; j < tam; ++j)
		{
			E = memmatr(A->n, A->m);
			cofacmmatrix(A, E, i, j, tam);
			signo2 = ((i + j) % 2 == 0) ? 1: -1;
			Adj->M[j][i] = (signo2) * (determmatrix(E, (tam - 1)));
		}
	freematr(E);
}

mmatrix* invermmatrix(mmatrix *A)
{
	mmatrix *Adj, *inverA;
	double determt;
	int i, j;

	if (A->n == A->m)
	{
		determt = determmatrix(A, A->n);
		inverA = memmatr(A->n, A->m);
		Adj = memmatr(A->n, A->m);
		adjmmatrix(A, Adj, A->n);
		for (i = 0; i < A->n; ++i)
			for (j = 0; j < A->n; ++j)
				inverA->M[i][j] = Adj->M[i][j]/determt;

		freematr(Adj);
	}
    return inverA;
}

double** inverse(double** matrix, int size){
	mmatrix *A = (mmatrix*)malloc(sizeof(mmatrix));
	A->M = matrix;
	A->n = size;
	A->m = size;
    int i;
	mmatrix *invA = invermmatrix(A);
    double** returnMatrix = (double**)calloc(size,sizeof(double*));
    for(i=0;i<size;i++)
        returnMatrix[i] = (double*)calloc(size,sizeof(double));

    copyAndSetPartOfArray2(returnMatrix, invA->M,size, size,0,0,0,0);
    freematr(invA);
    free(A);
	return returnMatrix;
}

double** multiply_MM(double** m, double** n, int r1,int c1, int r2, int c2)
{
    int i,j,k;

    double **result = (double**)calloc(r1,sizeof(double*));
    for(i=0;i<r1;i++)
        result[i] = (double*)calloc(c2,sizeof(double));

    for(i=0; i < r1; i++){
        for(j=0; j < c2; j++){
            for (k = 0; k < c1; k++)
                result[i][j] += m[i][k] * n[k][j];
        }
    }
    return result;
}

double** add_MM(double** m, double** n, int row, int col)
{
    int i,j;
    double **result = (double**)calloc(row,sizeof(double*));
    for(i=0;i<row;i++)
        result[i] = (double*)calloc(col,sizeof(double));

    for(i=0;i<row;i++)
        for(j=0;j<col;j++)
            result[i][j] = m[i][j] + n[i][j];

    return result;
}

double** create2dMatrix(int size)
{
     double** m = new double*[size];
     for(int i = 0; i < size; i++)
        m[i] = new double[size];

    return m;

}


void printNormalMatrix(double** m, int size)
{
    for(int i = 0; i < size; i++)
     {
         printf("\n");
         for(int j = 0; j < size; j++)
         {
            printf("%f - ", m[i][j]);
         }
     }
      printf("\n");
}

std::complex<double>** inverseComplexMatrix(std::complex<double>** matrix, int size)
{
    // first seperate real and imaginary matrices
    double** matrixReal = create2dMatrix(size);
    double** matrixImag = create2dMatrix(size);
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            matrixReal[i][j] = matrix[i][j].real();
            matrixImag[i][j] = matrix[i][j].imag();
        }
    }

    // copy real and imaginary matrices since we will be lost them when take inverse
    double** matrixRealCopy = create2dMatrix(size);
    double** matrixImagCopy = create2dMatrix(size);
    copyAndSetPartOfArray2(matrixRealCopy, matrixReal, size, size, 0, 0,  0, 0);
    copyAndSetPartOfArray2(matrixImagCopy, matrixImag, size, size, 0, 0,  0, 0);


    // myMatrix =  M = A + i*B so formulation is below
    // (A + i*B)^-1 = (A + B*A^-1*B)^-1 - i*(B + A*B^-1*A)^-1

    //A^-1
    double** invReal = inverse(matrixRealCopy, size);
    //B*A^-1
    double** bMulAinv = multiply_MM(matrixImag, invReal, size, size,  size,  size);
    //B*A^-1*B
    double** bMulAinvMulB = multiply_MM(bMulAinv, matrixImag,  size, size,  size,  size);
    //(A + B*A^-1*B)
    double** aAddbMulAinvMulB = add_MM(matrixReal, bMulAinvMulB, size, size);
    //(A + B*A^-1*B)^-1
    double** inverseOfReal = inverse(aAddbMulAinvMulB, size);

    delete[] invReal;
    delete[] bMulAinv;
    delete[] bMulAinvMulB;
    delete[] aAddbMulAinvMulB;


    //B^-1
    double** invImag = inverse(matrixImagCopy, size);
    //A*B^-1
    double** aMulBinv = multiply_MM(matrixReal, invImag,  size, size,  size,  size);
    //A*B^-1*A
    double** aMulBinvMulA = multiply_MM(aMulBinv, matrixReal,  size, size,  size,  size);
    //(B + A*B^-1*A)
    double** bAddaMulBinvMulA = add_MM(matrixImag, aMulBinvMulA, size, size);
    //(B + A*B^-1*A)^-1
    double** inverseOfImag = inverse(bAddaMulBinvMulA, size);

    delete[] invImag;
    delete[] aMulBinv;
    delete[] aMulBinvMulA;
    delete[] bAddaMulBinvMulA;

    delete[] matrixReal;
    delete[] matrixImag;
    delete[] matrixImagCopy;
    delete[] matrixRealCopy;


    std::complex<double>** result = new std::complex<double>*[size];
     for(int i = 0; i < size; i++)
     {
         result[i] = new std::complex<double>[size];
         for(int j = 0; j < size; j++)
         {
            result[i][j] = std::complex<double> (inverseOfReal[i][j],-1*inverseOfImag[i][j]);
         }
     }

    delete[] inverseOfReal;
    delete[] inverseOfImag;


    return result;

}