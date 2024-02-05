#include <iostream>
#include <ctime> 
#include <math.h>
#include <string>
#include <complex>
#include "./gnuplot-iostream.h"

using namespace std;

extern "C" void zheev_(char*, char*, int*, complex<double>*, int*, double*, complex<double>*, int*, double*, int*);
void calc_eigen(complex<double>* mat, int N_mat, double* eigen_out);

int main(){

    srand((unsigned)time(0)); // Seeding the random generator

    Gnuplot gp1;

    // We initilize arrays for our variables
    const int N = 500;
    complex<double> A[N][N];
    complex<double> H_sum[N][N];
    double eigenvals[N];
    double x_r, x_i;


    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){

            x_i = double(rand())/double(RAND_MAX);
            x_r = double(rand())/double(RAND_MAX);
            A[i][j] = complex<double> {x_i, x_r};
        }    
    }    

    // Adding the conjugate transpose manually (there are better ways but we ignore them for now)
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){

            H_sum[i][j] = A[i][j] + conj(A[j][i]);
        }    
    }  

    // Computing the eigenvalues and storing them in "eigenvals"
    calc_eigen(H_sum[0], N, eigenvals);

    // Ploting and exportnig
    gp1 << "set style fill solid border -1 \n";
    gp1 << "set key off \n";
    gp1 << "set title " << "\"A + A'  Eigenvalues \"" <<"\n";
    gp1 << "set key off "<<"\n";
    gp1 << "binwidth=1.5 "<<"\n";
    gp1 << "bin(x,width)=width*floor(x/width) "<<"\n";
    gp1 << "set output " << "\"hist-sum.png\"" << "\n"; // Comment this line and the next for plot only mode
    gp1 << "set terminal pngcairo size 2000, 2000 fontscale 4\n";
    gp1 << "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes \n";

    for (int i=0; i<N-1; i++){
        gp1 << to_string(eigenvals[i])<< "\n";
    }

    return 0;
}

void calc_eigen(complex<double>* mat, int N_mat, double* eigen_out){

    // Calculates eigenvalues of hermitian matrices
    // This functions uses the ZHEEV function from LAPACK lib.

    int LWORK = 6*N_mat;
    double RWORK[3*N_mat - 2];
    complex<double> WORK[6*N_mat];
    int info;

    char par1 = 'N';
    char par2 = 'U';

    zheev_(&par1, &par2, &N_mat, mat, &N_mat, eigen_out, WORK, &LWORK, RWORK, &info);

    cout<<"ZHEEV info: "<<info<<"\n";

}