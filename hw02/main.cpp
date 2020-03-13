// HW01 Gram-Schmidt Process 
// Name: Yu-Hsiu Huang     ID: 104061249
// Date: 20180306
#include <iostream>
#include <math.h>
#include "MAT.h"
#include "VEC.h"
using namespace std;

// function declaration
void LoadMatrix(MAT &A);                    // load data from .dat file
void Algo1(MAT &A, MAT &g);                 // perform algorithm 1.
void Algo2(MAT &A, MAT &g);                 // perform algorithm 2.
void Algo3(MAT &A, MAT &g);                 // perform algorithm 3.
void Algo4(MAT &A, MAT &g);                 // perform algorithm 4.
void Sigma(MAT &g);                         // calculate sigma

int main(void)
{
    int dim;                                // dimension of matrix

    cin >> dim;                             // load the matrix dimension
    MAT A(dim);                             // matrix A
    MAT g(dim);                             // transpose of G

    LoadMatrix(A);                          // load matrix A from .dat file

    // perform Gram-Schimdt Process: uncomment to perform the certain algoritm    
    Algo1(A, g);                            // perform algorithm 1.
    //Algo2(A, g);                            // perform algorithm 2.
    //Algo3(A, g);                            // perform algorithm 3.
    //Algo4(A, g);                            // perform algorithm 4.
    
    Sigma(g);                               // calculate Sigma
    return 0;
}
void LoadMatrix(MAT &A)                     // load data from .dat file
{
    int i, j;

    for (i=0; i<A.dim(); i++) {             // load A
        for (j=0; j<A.dim(); j++) {
            cin >> A[i][j];
        }
    }
}
void Algo1(MAT &A, MAT &g)                  // Algorithm 1
{
    int k, i;
    VEC zero(A.dim());                      // zero vector
    MAT a=(A.tpose());                      // transpose of A
    // Gram-Schmidt process starts
    g[0]=a[0];                              // g1 = a1
    for (k=1; k<A.dim(); k++) {
        g[k]=zero;
        for (i=0; i<k; i++) {
            g[k]+=((a[k]*g[i])*g[i]/(g[i]*g[i]));
        }
        g[k]=a[k]-g[k];
    }
    // process ends
    cout << "Algorithm 1. Gram-Schmidt Process" <<'\n';                        
}
void Algo2(MAT &A, MAT &g)                  // Algorithm 2
{
    int k, i;
    MAT a=(A.tpose());                      // transpose of A
    // Gram-Schmidt process starts
    g[0]=a[0];                              // g1 = a1
    for (k=1; k<A.dim(); k++) {
        g[k]=a[k];
        for (i=0; i<k; i++) {
            g[k]-=(g[k]*g[i]*g[i])/(g[i]*g[i]);
        }
    }
    // process ends
    cout << "Algorithm 2. Modified Gram-Schmidt Process, 1" <<'\n';
}
void Algo3(MAT &A, MAT &g)                  // Algorithm 3
{
    int k, i;
    MAT a=(A.tpose());                      // transpose of A
    // Gram-Schmidt process starts
    g[0]=a[0];                              // g1 = a1
    for (k=1; k<A.dim(); k++) {
        g[k]=a[k];
        for (i=0; i<k; i++) {
            g[k]-=(g[k]*g[i])/(g[i]*g[i])*g[i];
        }
    }
    // process ends
    cout << "Algorithm 3. Modified Gram-Schmidt Process, 2" <<'\n';
}
void Algo4(MAT &A, MAT &g)                  // Algorithm 4
{
    int i, j, k;
    double alpha;
    MAT a=(A.tpose());                      // transpose of A
    // Gram-Schmidt process starts
    g[0]=a[0];                              // g1 = a1
    for (i=0; i<A.dim(); i++) {
        g[i]=a[i];
    }
    for (j=0; j<(A.dim()-1); j++) {
        alpha=(g[j]*g[j]);                  // inner product of g[j]
        for (k=j+1; k<A.dim(); k++) {
            g[k]-=(g[k]*g[j]/alpha)*g[j];
        }
    }
    // process ends
    cout << "Algorithm 4. Modified Gram-Schmidt Process, 3" <<'\n';
}
void Sigma(MAT &g)                          // Sigma Calculation
{
    double s=0;
    int i, j;
    MAT G=(g.tpose());                      // G
    MAT D=(g*G);                            // D should be a diagonal matrix

    for (i=0; i<g.dim(); i++) {
        for (j=0; (j<g.dim()) and (j!=i); j++) {
            s+=(D[i][j]*D[i][j]);
        }
    }
    cout << "sigma = " << sqrt(s) << endl;  // output the square root
}
