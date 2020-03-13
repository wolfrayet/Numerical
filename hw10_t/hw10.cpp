//  EE4070 Numerical Analysis
//  HW 10. Nonlinear Network 
//  Yu-Hsiu Huang
#include <iostream>
#include <fstream>
#include <math.h>
#include "VEC.h"
#include "MAT.h"
#define R0 1
#define a 0.1
#define b 1
#define K 1
#define e 0.0000001
using namespace std;
void Q1stamp(MAT &A, VEC &F, VEC x, const int i, const int j) 
{   
                                    // Q1 stamp function
                                    // requirement: Vi > Vj for a resistor
    double v = x[i]-x[j];           // voltage across the resistor
    double r = R0 + a*fabs(v);      // voltage-dependent resistance
    F[i] -= v/r;
    F[j] += v/r;
    A[i][i] += R0/(r*r);
    A[i][j] -= R0/(r*r);
    A[j][i] -= R0/(r*r);
    A[j][j] += R0/(r*r);
}
void Q2stamp(MAT &A, VEC &F, VEC x, const int i, const int j, const int k)
{
                                    // Q2 stamp fuction for resistor k, where Vi Vj are the voltage at two ends
                                    // Vi > Vj
    double v = x[i]-x[j];           // voltage across the resistor
    double R = R0+K*x[k];           // temperature-dependent resistance
    F[i] -= v/R;
    F[j] += v/R;
    F[k] -= x[k] - b*v*v/R;
    A[i][i] += 1/R;
    A[i][j] -= 1/R;
    A[j][i] -= 1/R;
    A[j][j] += 1/R;
    A[i][k] += K*v/(R*R);
    A[j][k] -= K*v/(R*R);
    A[k][i] -= 2*b*v/R;
    A[k][j] += 2*b*v/R;
    A[k][k] += 1 + b*v*v*K/(R*R);
}
void Rstamp(MAT &A, VEC x, const int i, const int j) 
{   
                                    // stamp function for linear system
    A[i][i] += 1/R0;
    A[i][j] -= 1/R0;
    A[j][i] -= 1/R0;
    A[j][j] += 1/R0;
}
int main(void)
{
    fstream fout;
    fout.open("result.txt", fstream::out);  // the result is saved in result.txt
    double V;
/*  // the following code is for Q1
    double err;
    MAT A(9);                           // Jacobian Matrix
    VEC F(9); VEC x(9); VEC xt(9);
    fout << "Voltage" << " " << "It" << " " << "I2" << " " << "I7" << " " << "I12" << endl;
    for (V=0; V<=5; V+=0.1) {           // scan from 0V to 5V, increment=0.1V
        err=1+e;
        while (err>=e) {                // iteration start
            // initislize the linear system
            for (int i=0; i<9; i++) {
                F[i]=0; xt[i]=0;
                for (int j=0; j<9; j++) A[i][j]=0;
            }
            // construct the linear system for Newton Method
            Q1stamp(A, F, x, 0, 1);
            Q1stamp(A, F, x, 0, 3);
            Q1stamp(A, F, x, 1, 2);
            Q1stamp(A, F, x, 1, 4);
            Q1stamp(A, F, x, 2, 5);
            Q1stamp(A, F, x, 3, 4);
            Q1stamp(A, F, x, 3, 6);
            Q1stamp(A, F, x, 4, 5);
            Q1stamp(A, F, x, 6, 7);
            Q1stamp(A, F, x, 7, 4);
            Q1stamp(A, F, x, 7, 8);
            Q1stamp(A, F, x, 8, 5);
            for (int j=0; j<9; j++) {   // set the voltage source node and ground node
                if (j==0) A[0][j]=1;
                else A[0][j]=0;
                if (j==5) A[5][j]=1;
                else A[5][j]=0;
            }
            F[0]=-x[0]+V; F[5]=-x[5];
            // LU decomposition to implement the Newton method
            A=luFact(A);
            xt=fwdSubs(A, F);
            xt=bckSubs(A, xt);
            x+=xt;
            // error estimate
            err=norms(F, 2);
        }
        // output the current through resistors r2, r7, and r12 
        fout << V << " ";
        fout << (x[0]-x[1])/(R0+a*fabs(x[0]-x[1])) + (x[0]-x[3])/(R0+a*fabs(x[0]-x[3])) << " ";
        fout << (x[1]-x[2])/(R0+a*fabs(x[1]-x[2]))  << " ";
        fout << (x[4])/(R0+a*fabs(x[4])) << " ";
        fout << (x[7]-x[8])/(R0+a*fabs(x[7]-x[8])) << endl;
    }
*/

    // the following code is for Q2
    double err;
    MAT A(21);                          // Jacobian Matrix
    VEC F(21); VEC x(21); VEC xt(21);
    fout << "Voltage" << " " << "T2" << " " << "T7" << " " << "T12" << endl;
    for (V=0; V<=5; V+=0.1) {           // scan from 0V to 5V, increment=0.1V
        err=1+e;
        while (err>=e) {                // iteration start
            // initialize the linear system
            for (int i=0; i<21; i++) {
                F[i]=0; xt[i]=0;
                for (int j=0; j<21; j++) A[i][j]=0;
            }
            // construct the linear system for Newton Method
            Q2stamp(A, F, x, 0, 1, 9);
            Q2stamp(A, F, x, 1, 2, 10);
            Q2stamp(A, F, x, 0, 3, 11);
            Q2stamp(A, F, x, 1, 4, 12);
            Q2stamp(A, F, x, 2, 5, 13);
            Q2stamp(A, F, x, 3, 4, 14);
            Q2stamp(A, F, x, 4, 5, 15);
            Q2stamp(A, F, x, 3, 6, 16);
            Q2stamp(A, F, x, 7, 4, 17);
            Q2stamp(A, F, x, 8, 5, 18);
            Q2stamp(A, F, x, 6, 7, 19);
            Q2stamp(A, F, x, 7, 8, 20);
            for (int j=0; j<21; j++) {  // set the voltage source node and ground node
                if (j==0) A[0][j]=1;
                else A[0][j]=0;
                if (j==5) A[5][j]=1;
                else A[5][j]=0;
            }
            F[0]=-x[0]+V; F[5]=-x[5];
            // LU decomposition to implement the Newton method
            A=luFact(A);
            xt=fwdSubs(A, F);
            xt=bckSubs(A, xt);
            x+=xt;
            // error estimate
            err=norms(F, 2);
        }
        // output the temperature of the resistor r2, r7, r12
        fout << V << " " << x[10] << " " << x[15] << " " << x[20] << endl;
    }

/*
    // the following code is for the linear resistor network of the same dimension
    MAT A(9);
    VEC x(9); VEC y(9); VEC F(9);
    fout << "Voltage" << " " << "It" << " " << "I2" << " " << "I7" << " " << "I12" << endl;
    for (V=0; V<5; V+=0.1) {
        for (int i=0; i<9; i++) {
                F[i]=0;
                for (int j=0; j<9; j++) A[i][j]=0;
        }
        Rstamp(A, x, 0, 1);
        Rstamp(A, x, 0, 3);
        Rstamp(A, x, 1, 2);
        Rstamp(A, x, 1, 4);
        Rstamp(A, x, 2, 5);
        Rstamp(A, x, 3, 4);
        Rstamp(A, x, 3, 6);
        Rstamp(A, x, 4, 5);
        Rstamp(A, x, 6, 7);
        Rstamp(A, x, 7, 4);
        Rstamp(A, x, 7, 8);
        Rstamp(A, x, 8, 5);
        for (int j=0; j<9; j++) {   // set the voltage source node and ground node
            if (j==0) A[0][j]=1;
            else A[0][j]=0;
            if (j==5) A[5][j]=1;
            else A[5][j]=0;
        }
        F[0]=V;
        // LU decomposition 
        A=luFact(A);
        y=fwdSubs(A, F);
        x=bckSubs(A, y);
        fout << V << " " << (x[0]-x[1])/R0 + (x[0]-x[3])/R0;
        fout << " " << (x[1]-x[2])/R0;
        fout << " " << (x[4]-x[5])/R0;
        fout << " " << (x[7]-x[8])/R0 <<endl;
    }
*/
    fout.close();
    return 0;
}
