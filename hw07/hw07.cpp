//  EE4070 Numerical Analysis
//  HW7. Matrix Eigenvalues
//  Name: Yu-Hsiu Huang     ID: 104061249
//  Date: 2018/04/26

#include <iostream>
#include "VEC.h"
#include "MAT.h"
#define MAXITER 100000
#define TOL 0.000000001 
#define MU 0.5
using namespace std;

void Answer(MAT A);
int main(void)
{
    int dim=0;
    int step=0;
    cin >> dim;         // read the dimension of the matrix
    MAT A(dim);         // construct a matrix
    // == read the matrix from data file ==
    for (int i=0; i<dim; i++) {         
        for (int j=0; j<dim; j++) {
            cin >> A[i][j];
        }
    }
    // ======== find eigenvalues =========
    //step=EVqr(A, TOL, MAXITER);
    step=EVqrshifted(A, MU, TOL, MAXITER);
    // ========== print answer ===========
    cout << "Dimension of matrix: " << dim <<'\n';
    cout << "Step: " << step <<'\n';
    Answer(A);
    return 0;
}
void Answer(MAT A)
{
    int dim=A.dim();
    double *EV=(double *)calloc(A.dim(), sizeof(double));
    for (int i=0; i<A.dim(); i++) EV[i]=A[i][i];
    sort <double*> (EV, EV+dim);
    cout << "Three Largest Eigenvalues: " << '\n';
    for (int i=A.dim()-1; i>A.dim()-4; i--) cout << EV[i] << " ";
    cout << '\n' << "Three Smallest Eigenvalues: " << '\n';
    for (int i=0; i<3; i++) cout << EV[i] << " ";
    cout << '\n';
    free(EV);
}