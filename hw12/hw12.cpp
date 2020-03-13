//  EE4070 Numerical Analysis
//  HW 12. RLC Circuit II
//  Yu-Hsiu Huang   104061249
//  3 June, 2018
#include <iostream>
#include <fstream>
#include "VEC.h"
#include "MAT.h"
using namespace std;
#define R   1
#define Rs  0.1
#define C   1
#define Cs  0.1
#define R   1
#define V   1
#define L   1
#define h   0.2
void TrapeSet(MAT &M)                   // Set matrix for Trapezoidal Method
{
    M[0][0] = 1 + h/2/Rs/Cs;    M[0][1] = 0;            M[0][2] = 0;        M[0][3] = h/2/Cs;
    M[1][0] = -1;               M[1][1] = 1;            M[1][2] = 0;        M[1][3] = R;
    M[2][0] = 0;                M[2][1] = 0;            M[2][2] = 1;        M[2][3] = -h/2/C;
    M[3][0] = 0;                M[3][1] = -h/2/L;       M[3][2] = h/2/L;    M[3][3] = 1;
    M = luFact(M);                      // decompose the matrix into LU 
}
void Trape(VEC &x, MAT &M)              // Trapezoidal Method
{
    VEC xt(4);
    xt[0] = x[0]*(1 - h/2/Rs/Cs) - x[3]*h/2/Cs + V*h/Rs/Cs; 
    xt[1] = 0; 
    xt[2] = x[2] + x[3]*h/2/C;
    xt[3] = x[1]*h/2/L - x[2]*h/2/L + x[3];
    xt = fwdSubs(M, xt);                // solve x via LU decomposition
    x = bckSubs(M, xt);
}
void GearSet(MAT &M)                    // Set matrix for 2nd order Gear Method 
{
    M[0][0] = 1 + 2*h/3/Rs/Cs;  M[0][1] = 0;            M[0][2] = 0;        M[0][3] = 2*h/3/Cs;
    M[1][0] = -1;               M[1][1] = 1;            M[1][2] = 0;        M[1][3] = R;
    M[2][0] = 0;                M[2][1] = 0;            M[2][2] = 1;        M[2][3] = -2*h/3/C;
    M[3][0] = 0;                M[3][1] = -2*h/3/L;     M[3][2] = 2*h/3/L;  M[3][3] = 1;
    M = luFact(M);                      // decompose the matrix into LU
}
void Gear(VEC &x, VEC &xtt, MAT &M)     // 2nd order Gear Method
{
    VEC xt(4);
    xt[0] = x[0]*4/3 - xtt[0]/3 + 2*V*h/3/Rs/Cs; 
    xt[1] = 0; 
    xt[2] = x[2]*4/3 - xtt[2]/3;
    xt[3] = x[3]*4/3 - xtt[3]/3;
    xtt = x;                            // save the t-h terms
    xt = fwdSubs(M, xt);                // solve x via LU decomposition
    x = bckSubs(M, xt);
}
int main(void)
{
    fstream fout;
    fout.open("result.txt", fstream::out);  // create a txt file to record the result
    fout << "time" << " " << "V1" << " " << "V2" << " "  << "V3" << " " << "I" <<endl;
    VEC x(4), max(4), min(4);               // vector for values, maximum, minimum
    VEC xtt(4);                             // the temparol terms for 2nd Gear Method
    MAT M(4);                               // matrix for system
    x[0]=0; x[1]=0; x[2]=0; x[3]=0;         // initial values
    max = x; min = x; xtt = x;
/*
    fout << "0" << " " << x[0] << " " << x[1] << " " << x[2]  << " " << x[3] <<endl;
    // ========= Trapezoidal =========
    TrapeSet(M);                            // set the matrix for the system
    for (double t=h; t<=10; t+=h) {         // solve the differential system from t=0 to t=10
        Trape(x, M);
        // === write the result to the txt file ===
        fout << t << " " << x[0] << " " << x[1] << " " << x[2]  << " " << x[3] <<endl;
        // ====== find the maximum and minimum ====
        if (x[0] > max[0]) max[0] = x[0];
        else if (x[0] < min[0]) min[0] = x[0];
        if (x[1] > max[1]) max[1] = x[1];
        else if (x[1] < min[1]) min[1] = x[1];
        if (x[2] > max[2]) max[2] = x[2];
        else if (x[2] < min[2]) min[2] = x[2];
        if (x[3] > max[3]) max[3] = x[3];
        else if (x[3] < min[3]) min[3] = x[3];
    }
*/

/**/
    // =========== Gear 2nd order ==========
    fout << "0" << " " << x[0] << " " << x[1] << " " << x[2]  << " " << x[3] <<endl;
    TrapeSet(M);                            // use trapezoidal method to find x(h)
    Trape(x, M);
    fout << h << " " << x[0] << " " << x[1] << " " << x[2]  << " " << x[3] <<endl;
    if (x[0] > max[0]) max[0] = x[0];
    else if (x[0] < min[0]) min[0] = x[0];
    if (x[1] > max[1]) max[1] = x[1];
    else if (x[1] < min[1]) min[1] = x[1];
    if (x[2] > max[2]) max[2] = x[2];
    else if (x[2] < min[2]) min[2] = x[2];
    if (x[3] > max[3]) max[3] = x[3];
    else if (x[3] < min[3]) min[3] = x[3];
    GearSet(M);                             // set the matrix for the system
    for (double t=2*h; t<=10; t+=h) {
        // === write the result to the txt file ===
        Gear(x, xtt, M);                    // solve the differential system from t=0 to t=10
        fout << t << " " << x[0] << " " << x[1] << " " << x[2]  << " " << x[3] <<endl;
        // ====== find the maximum and minimum ====
        if (x[0] > max[0]) max[0] = x[0];
        else if (x[0] < min[0]) min[0] = x[0];
        if (x[1] > max[1]) max[1] = x[1];
        else if (x[1] < min[1]) min[1] = x[1];
        if (x[2] > max[2]) max[2] = x[2];
        else if (x[2] < min[2]) min[2] = x[2];
        if (x[3] > max[3]) max[3] = x[3];
        else if (x[3] < min[3]) min[3] = x[3];
    }
/**/
    fout.close();
    // ============= output ======================
    cout << "max V1 V2 V3 I" <<endl;
    for (int i=0; i<4; i++) cout << " " << max[i];
    cout << endl;
    cout << "min V1 V2 V3 I" <<endl;
    for (int i=0; i<4; i++) cout << " " << min[i];
    cout << endl;
    return 0;
}
