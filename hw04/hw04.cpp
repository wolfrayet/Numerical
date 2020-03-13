/*  EE4070 Numerical Analysis
    HW4. Linear Iterative Methods
    Name: Yu-Hsiu Huang     ID: 104061249
    Date: 2018/03/30
*/
#include <iostream>
#include "VEC.h"
#include "MAT.h"
#include <math.h>
using namespace std;

void Construct(MAT &A, VEC &b, double g);
void Symmetrize(MAT &A, VEC &b, double g);
void GetAnswer(VEC x, int const step, int const maxIter);

int main(void) 
{
    int dim=441;                            // number of total nodes
    int step=0;                             // iterative steps
    int maxIter=16000;
    double g=0.01;                          // conductance
    double tol=0.0000001;                   // tolerance
    MAT A(dim);
    VEC b(dim);
    VEC x(dim);
    x[10]=1;
    x[230]=0;
    Construct(A, b, g);                     // construct the linear system
    Symmetrize(A, b, g);                    // symmetrize linear system
    //step=jacobi(A, b, x, maxIter, tol);
    //step=gaussSeidel(A, b, x, maxIter, tol);
    step=sgs(A, b, x, maxIter, tol);
    GetAnswer(x, step, maxIter);
    return 0;
}
void Construct(MAT &A, VEC &b, double g)    // construct the linear system
{
    for (int i=0; i<b.len(); i++) {
        if (i!=10) b[i]=0;
        else b[i]=1;
    }
    for (int i=0; i<20; i++) {
        for (int j=0; j<=20; j++) {
            int n=i+j*21;
            if (i==10 && j==0) {
                A[n][n]=1;
                A[n+1][n]-=g;
                A[n+1][n+1]+=g;
            }
            else if (i+1==10 && j==0) {
                A[n+1][n+1]=1;
                A[n][n+1]-=g;
                A[n][n]+=g;
            }
            else if (i==19 && j==10) {
                A[n+1][n+1]=1;
                A[n][n+1]-=g;
                A[n][n]+=g;
            }
            else {
                A[n+1][n+1]+=g;
                A[n][n]+=g;
                A[n][n+1]-=g;
                A[n+1][n]-=g;
            }
        }
    }
    for (int i=0; i<b.len()-21; i++) {
        if (i==10) {
            A[i][i]=1;
            A[i+21][i]-=g;
            A[i+21][i+21]+=g;
        }
        else if (i==230) {
            A[i][i]=1;
            A[i+21][i]-=g;
            A[i+21][i+21]+=g;
        }
        else if (i+21==230) {
            A[i+21][i+21]=1;
            A[i][i+21]-=g;
            A[i][i]+=g;
        }
        else {
            A[i][i]+=g;
            A[i+21][i+21]+=g;
            A[i][i+21]-=g;
            A[i+21][i]-=g;
        }
    }

}
void Symmetrize(MAT &A, VEC &b, double g)   // symmetrize linear system
{
    for (int i=0; i<A.dim(); i++) {
        if (i!=10 && A[i][10]!=0) {
            A[i][10]=0;
            b[i]+=g;
        }
        else if (i!=230 && A[i][230]!=0) {
            A[i][230]=0;
        }
        else;
    }
}
void GetAnswer(VEC x, int const step, int const maxIter)
{
    if (step<maxIter) {
        cout << "The iteration steps: " << step << '\n';
        cout << "Voltage value for the south-west corner: " << x[20] << '\n';
        cout << "Voltage value for the north-east corner: " << x[420] << '\n';
        cout << "Voltage value for the center of east side: " << x[430] << '\n';
    }
    else {
        cout << "Not convergent." << '\n';
    }
}