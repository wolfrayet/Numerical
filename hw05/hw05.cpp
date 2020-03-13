/*  EE4070 Numerical Analysis
    HW05. Conjugate Gradient Methods
    Name: Yu-Hsiu Huang  ID: 104061249
    Date: 2018/04/04
*/
#include <iostream>
#include <iomanip>
#include "VEC.h"
#include "MAT.h"
using namespace std;

void Construct(MAT &A, VEC &b, int const num, double g, int source[2], int ground[2]);
void Symmetrize(MAT &A, VEC &b, int const num, double g, int source[2], int ground[2]);
void LU(MAT &A, VEC b, VEC &x);
void Answer(VEC x, int source[2], double g, int step, int num, int maxIter);

int main(void)
{
    int num=100;
    int dim=(num+1)*(num+1);
    int step=0;
    int maxIter=dim;
    int source[2];                      // location of source, source[0] for verticle, source[0] for horizon
    int ground[2];                      // location of ground, ground[0] for verticle, ground[0] for horizon
    double g=(double)num/2000;
    double tol=0.000000000001;

    source[0]=num/2;                    // the location of node connected to the source
    source[1]=0;
    ground[0]=num;                      // the location of grounded node
    ground[1]=num/2;
    
    MAT A(dim);
    VEC b(dim);
    VEC x(dim);
    
    Construct(A, b, num, g, source, ground);        // Construct the linear system
    Symmetrize(A, b, num, g, source, ground);       // Symmetrize the linear system

    //step=cg(A, b, x, maxIter, tol);                 // Conjugate Gradient Method
    step=cg2(A, b, x, maxIter, tol);                // Conjugate Gradient Method, 2nd form
    //step=gaussSeidel(A, b, x, maxIter, tol);        // Gauss-Seidel Method
    //LU(A, b, x);                                    // LU decomposition

    Answer(x, source, g, step, num, maxIter);       // compute and output answers
    return 0;
}
void Construct(MAT &A, VEC &b, int const num, double g, int source[2], int ground[2])    // construct and initialize the linear system
{
    // the right-hand side of the linear system
    for (int i=0; i < b.len(); i++) {     
        if (i!=source[0]+source[1]*(num+1)) b[i]=0;
        else b[i]=1;
    }
    // the verticle resistors
    for (int i=0; i < num; i++) {
        for (int j=0; j <= num; j++) {
            int n=i+j*(num+1);
            if (i==source[0] && j==source[1]) {
                A[n][n]=1;
                A[n+1][n]-=g;
                A[n+1][n+1]+=g;
            }
            else if (i+1==source[0] && j==source[1]) {
                A[n+1][n+1]=1;
                A[n][n+1]-=g;
                A[n][n]+=g;
            }
            else if (i==ground[0] && j==ground[1]) {
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
    // the horizontal resistors
    for (int i=0; i < b.len()-num-1; i++) {
        if (i==source[0]+source[1]*(num+1)) {
            A[i][i]=1;
            A[i+num+1][i]-=g;
            A[i+num+1][i+num+1]+=g; 
        }
        else if (i+num+1==source[0]+source[1]*(num+1)) {
            A[i+num+1][i+num+1]=1;
            A[i][i+num+1]-=g;
            A[i][i]+=g;
        }
        else if (i==ground[0]+ground[1]*(num+1)) {
            A[i][i]=1;
            A[i+num+1][i]-=g;
            A[i+num+1][i+num+1]+=g;
        }
        else if (i+num+1==ground[0]+ground[1]*(num+1)) {
            A[i+num+1][i+num+1]=1;
            A[i][i+num+1]-=g;
            A[i][i]+=g;
        }
        else {
            A[i][i]+=g;
            A[i+num+1][i+num+1]+=g;
            A[i][i+num+1]-=g;
            A[i+num+1][i]-=g;
        }
    }
}
void Symmetrize(MAT &A, VEC &b, int const num, double g, int source[2], int ground[2])   // symmetrize linear system
{
    int node_s=source[0]+source[1]*(num+1);     // node connected to source
    int node_g=ground[0]+ground[1]*(num+1);     // node connected to ground

    for (int i=0; i < A.dim(); i++) {
        if (i!=node_s && A[i][node_s]!=0) {
            A[i][node_s]=0;
            b[i]+=g;
        }
        else if (i!=node_g && A[i][node_g]!=0) {
            A[i][node_g]=0;
        }
        else;
    }
}
void LU(MAT &A, VEC b, VEC &x)                  // LU decomposition
{
    VEC y(x.len());

    A=luFact(A);
    y=fwdSubs(A, b);                            // forward substitution
    x=bckSubs(A, y);                            // backward substitution
}
void Answer(VEC x, int source[2], double g, int step, int num, int maxIter) // output answer
{
    double current=0;
    int node_s=source[0]+source[1]*(num+1);     // the node connected to the source
    int side=0;                                 // source coneected to the east side or west side

    if (source[1]!=0) side=-num-1;
    else side=num+1;
    current=g*(3*x[node_s]-x[node_s-1]-x[node_s+1]-x[node_s+side]);
    cout << "The Solution for the network with number of side-resistor: " << num << '\n';
    if (step < maxIter) {
        cout << "The iteration steps: " << step << '\n';
        cout << "Voltage value for the south-west corner: " << setprecision(9)<<x[num] << '\n';
        cout << "Voltage value for the north-east corner: " << x[x.len()-num-1] << '\n';
        cout << "Voltage value for the center of east side: " << x[x.len()-num/2-1] << '\n';
        cout << "Equivalent Resistance: " << setprecision(11)<< 1/current << '\n';
    }
    else {
        cout << "The iterative method is not convergent." << '\n';
    }
}
