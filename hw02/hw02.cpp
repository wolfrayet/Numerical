/*  HW02. LU Decomposition
    Name: Yu-Hsiu Huang ID: 104061249
    Date: 2018/03/16
*/
#include <iostream>
#include "MAT.h"
#include "VEC.h"
using namespace std;

void Load(MAT &A, VEC &b);
void Print(VEC &b);
void Check(MAT &O, VEC &x, VEC &b);

int main(void)
{
    int dim;                    // dimension of matrix    
    cin >> dim;                 // load the matrix dimension
    MAT A(dim);                 // nonsingular matrix A
    MAT O(dim);                 // the original matrix
    VEC b(dim);                 // right-hand side vector b
    VEC y(dim);                 // save temporal solution to L
    VEC x(dim);                 // unknown vector x
    Load(A, b);                 // load the linear system parameters
    O=A;                        // save the original A
    A=luFact(A);                // in-place LU decomposition
    y=fwdSubs(A, b);            // forward substitution: Ly=b
    x=bckSubs(A, y);            // backward substitution: Ux=y
    Check(O, x, b);             // check whether x is root of the input linear system
    //Print(x);                   // uncomment to print vector x
    return 0;
}
void Load(MAT &A, VEC &b)       // load the parameters of linear system
{      
    for (int i=0; i<A.dim(); i++) {     // load matrix A
        for (int j=0; j<A.dim(); j++) {
            cin >> A[i][j];
        }
    }
    for (int i=0; i<b.len(); i++ ) {    // load right-hand side vector b
        cin >> b[i];
    }
 }
 void Print(VEC &b)             // output vector
 {
     for (int i=0; i<b.len(); i++) {
         cout << b[i] << '\n';
     }
 }
 void Check(MAT &A, VEC &x, VEC &b)     // check whether x is the root
 {
     VEC t(x.len());            // temparol vector 
     int check;                 // correctness 
     double error;              // error 
     t=A*x;
     check=0;                   // number of incorrect element
     // compare the erro between t and b vector, the error cannot be neglected if error > 0.01%
     for (int i=0; i<x.len(); i++) {
         error=(t[i]-b[i])/t[i];
         if (error>0.0001 || error<-0.0001) {
             cout << i+1 << "-th element of x is wrong." << '\n';
             check++;
         }
     }
     // output whether the root is correct or not
     if (check==0) {
         cout << "x is correct." << '\n';
     }
     else {
         cout << "x is wrong." << '\n';
     }
 }
