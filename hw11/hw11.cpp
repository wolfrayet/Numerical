//  EE4070 Numerical Analysis
//  HW11 RLC Circuit
//  Yu-Hsiu Huang   104061249
//  26 May, 2018
#include <iostream>
#include <fstream>
#include <math.h>
#include "VEC.h"
#include "MAT.h"
#define R 1
#define L 1
#define C 1
#define V 1
#define h 0.01
using namespace std;

void fwdEuler(VEC &x)       // forwrd Euler method for one step
{
    VEC xt(3);
    for (int i=0; i<3; i++) xt[i] = x[i];
    x[0] = V - xt[2]*R;
    x[1] = xt[1] + xt[2]*h/C;
    x[2] = (xt[0]-xt[1])*h/L + xt[2];
}
void bckEulerSet(MAT &M)    // setup the matrix for backward Euler method
{
    M[0][0] = h/L; M[0][1] = -h/L; M[0][2] = -1;
    M[1][0] = 0; M[1][1] = 1; M[1][2] = -h/C;
    M[2][0] = 1; M[2][1] = 0; M[2][2] = R;
    M = luFact(M);
}
void TrapeSet(MAT &M)
{
    M[0][0] = h/(2*L); M[0][1] = -h/(2*L); M[0][2] = -1;
    M[1][0] = 0; M[1][1] = 1; M[1][2] = -h/(2*C);
    M[2][0] = 1; M[2][1] = 0; M[2][2] = R;
    M = luFact(M);
}
void bckEuler(VEC &x, MAT &M)   // backward Euler method
{
    VEC xt(3);
    xt[0] = -x[2]; xt[1] = x[1]; xt[2] = V;
    xt = fwdSubs(M, xt);
    x = bckSubs(M, xt);
}
void Trape(VEC &x, MAT &M)
{
    VEC xt(3);
    xt[0] = -x[0]*h/(2*L) + x[1]*h/(2*L) - x[2]; 
    xt[1] = x[1] + x[2]*h/(2*C); 
    xt[2] = V;
    xt = fwdSubs(M, xt);
    x = bckSubs(M, xt);
}
int main(void)
{
    fstream fout;
    fout.open("result.txt", fstream::out);  // create a txt file to record the result
    fout << "time" << " " << "V1" << " " << "V2" << " " << "I" <<endl;
    VEC x(3), max(3), min(3);
    MAT M(3);
    x[0]=1; x[1]=0; x[2]=0;
    // ========= Setup Matrix =========
    //bckEulerSet(M);         // setup the matrix for the backward euler method
    TrapeSet(M);
    // ========== Evaluation ==========
    for (double t=0; t<=10; t+=h) {
        if (t==h) {
            max[0]=x[0]; max[1]=x[1]; max[2]=x[2];
            min[0]=x[0]; min[1]=x[1]; min[2]=x[2];
        }
        //fwdEuler(x);        // Forward Euler Method
        //bckEuler(x, M);     // Backword Euler Method
        Trape(x, M);
        // ====== find the maximum and minimum ====
        if (x[0] > max[0]) max[0] = x[0];
        else if (x[0] < min[0]) min[0] = x[0];
        if (x[1] > max[1]) max[1] = x[1];
        else if (x[1] < min[1]) min[1] = x[1];
        if (x[2] > max[2]) max[2] = x[2];
        else if (x[2] < min[2]) min[2] = x[2];

        // === write the result to the txt file ===
        fout << t << " " << x[0] << " " << x[1] << " " << x[2] <<endl;
    }
    fout.close();
    // ========= Print Result =========
    cout << "max V1 V2 I" <<endl;
    for (int i=0; i<3; i++) cout << " " << max[i];
    cout << endl;
    cout << "min V1 V2 I" <<endl;
    for (int i=0; i<3; i++) cout << " " << min[i];
    cout << endl;
    
    return 0;
}
