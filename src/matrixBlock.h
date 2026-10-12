#include <complex>

using namespace std;

class matrixBlock{
    public:
    int sizerow;
    int sizecol;
    int* indrows;
    int* indcols;
    complex<double>** submatrix;

    //class constructor with null elements
    matrixBlock(){} 
        
    //initialization of sizerow and allocation of submatrix
    void initmatrixBlock(int matsizerow, int matsizecol){
        sizerow=matsizerow;
        sizecol=matsizecol;
        indrows = new int[sizerow];
        indcols = new int[sizecol];
        submatrix = new complex<double>*[sizerow];
        for(int i=0; i<sizerow; i++){
            submatrix[i]=new complex<double>[sizecol];
        }
    }

    //initialization of submatrix
    void initsubmatrix(complex<double> initval){
        for(int i=0;i<sizerow; i++){
            for(int j=0;j<sizecol; j++){
                submatrix[i][j]=initval;
            }
        }
    }

    //assign value to submatrix
    void assignvalue(int indcol, int indrow, complex<double> value){
        submatrix[indcol][indrow]=value;
    }

    //class destructor
    ~matrixBlock(){
        //printf("Deleting matrixBlock\n"); 
        delete[] indrows;
        delete[] indcols;
        for(int i=0; i<sizerow; i++){
            delete[] submatrix[i];
        }
        delete[] submatrix;
    }

};