//  EE4070 Numerical Analysis
//  HW10 Nonlinear Resistor Networks
//  Name: Yu-Hsiu Huang ID: 104061248
//  Date: 2018/05/20
#include <iostream>
#include <fstream>
#include <math.h>
#include "VEC.h"
#include "MAT.h"
#define R 1
#define k 1
#define a 0.1
#define b 1
#define gnd 0
#define e 0.0000001
using namespace std;

double I1(double v1, double v2)
{
    double v=v1-v2;
    return v/(R + a*fabs(v));
}
void F1(VEC &x, VEC &F, double V)   // trnasformation function for Q1
{
    F[0] = I1(V, x[0]) + I1(x[1], x[0]) + I1(x[3], x[0]);
    F[1] = I1(x[0], x[1]) + I1(gnd, x[1]);
    F[2] = I1(V, x[2]) + I1(x[3], x[2]) + I1(x[4], x[2]);
    F[3] = I1(x[0], x[3]) + I1(x[2], x[3]) + I1(gnd, x[3]) + I1(x[5], x[3]);
    F[4] = I1(x[2], x[4]) + I1(x[5], x[4]);
    F[5] = I1(x[3], x[5]) + I1(x[4], x[5]) + I1(x[6], x[5]);
    F[6] = I1(gnd, x[6]) + I1(x[5], x[6]);
}
double diff(double v1, double v2, bool before) // differential for Q1
{
    double x=R+a*fabs(v1-v2);
    if (before == 1) {
        if (v1>v2) return (x-a*(v1-v2))/(x*x);
        else return (x+a*(v1-v2))/(x*x);
    }
    else {
        if (v1>v2) return (-x+a*(v1-v2))/(x*x);
        else return (-x-a*(v1-v2))/(x*x);
    }
}
void J1(VEC &x, MAT &J, double V)   // Jacobian for Q1
{
    J[0][0]=diff(V, x[0], 0) + diff(x[1], x[0], 0) + diff(x[3], x[0], 0);
    J[0][1]=diff(x[1], x[0], 1);
    J[0][3]=diff(x[1], x[0], 1);
    J[1][0]=diff(x[0], x[1], 1);
    J[1][1]=diff(x[0], x[1], 0) + diff(gnd, x[1], 0);
    J[2][2]=diff(V, x[2], 0) + diff(x[4], x[2], 0) + diff(x[3], x[2], 0);
    J[2][3]=diff(x[3], x[2], 1);
    J[2][4]=diff(x[4], x[2], 1);
    J[3][0]=diff(x[0], x[3], 1);
    J[3][2]=diff(x[2], x[3], 1);
    J[3][3]=diff(x[0], x[3], 0) + diff(x[2], x[3], 0) + diff(x[5], x[3], 0) + diff(gnd, x[3], 0);
    J[3][5]=diff(x[5], x[3], 1);
    J[4][2]=diff(x[2], x[4], 1);
    J[4][4]=diff(x[2], x[4], 0) + diff(x[5], x[4], 0);
    J[4][5]=diff(x[5], x[4], 1);
    J[5][3]=diff(x[3], x[5], 1);
    J[5][4]=diff(x[4], x[5], 1);
    J[5][5]=diff(x[3], x[5], 0) + diff(x[4], x[5], 0) + diff(x[6], x[5], 0);
    J[5][6]=diff(x[6], x[5], 1);
    J[6][5]=diff(x[5], x[6], 1);
    J[6][6]=diff(x[5], x[6], 0) + diff(gnd, x[6], 0);
}
double I(double v1, double v2, double T)    // current for Q2
{
    return (v1-v2)/(R+k*T);
}
double Ti(double v1, double v2, double T)   // temperature function for Q2
{
    return b*(v1-v2)*(v1-v2)/(R+k*T);
}
void F2(VEC &x, VEC &F, double V)   // transformation function for Q2
{
    F[0]=I(V, x[0], x[7]) + I(x[1], x[0], x[8]) + I(x[3], x[0], x[10]);
    F[1]=I(x[0], x[1], x[8]) + I(gnd, x[1], x[11]);
    F[2]=I(V, x[2], x[9]) + I(x[3], x[2], x[12]) + I(x[4], x[2], x[14]);
    F[3]=I(x[0], x[3], x[10]) + I(x[2], x[3], x[12]) + I(gnd, x[3], x[13]) + I(x[5], x[3], x[15]);
    F[4]=I(x[2], x[4], x[14]) + I(x[5], x[4], x[17]);
    F[5]=I(x[3], x[5], x[15]) + I(x[4], x[5], x[17]) + I(x[6], x[5], x[18]);
    F[6]=I(gnd, x[6], x[16]) + I(x[5], x[6], x[18]);
    F[7]=x[7]-Ti(V, x[0], x[7]);
    F[8]=x[8]-Ti(x[0], x[1], x[8]);
    F[9]=x[9]-Ti(V, x[2], x[9]);
    F[10]=x[10]-Ti(x[0], x[3], x[10]);
    F[11]=x[11]-Ti(x[1], gnd, x[11]);
    F[12]=x[12]-Ti(x[2], x[3], x[12]);
    F[13]=x[13]-Ti(x[3], gnd, x[13]);
    F[14]=x[14]-Ti(x[2], x[4], x[14]);
    F[15]=x[15]-Ti(x[3], x[5], x[15]);
    F[16]=x[16]-Ti(gnd, x[6], x[16]);
    F[17]=x[17]-Ti(x[4], x[5], x[17]);
    F[18]=x[18]-Ti(x[5], x[6], x[18]);
}
void stamp(MAT &J, VEC &x, int i, int j, int t) // resistance and temperature stamping for Q2 
{
    J[i][i]-=1/(R+k*x[t]);
    J[i][j]+=1/(R+k*x[t]);
    J[j][i]+=1/(R+k*x[t]);
    J[j][j]-=1/(R+k*x[t]);
    J[i][t]=(x[i]-x[j])*k/((R+k*x[t])*(R+k*x[t]));
    J[j][t]=-(x[i]-x[j])*k/((R+k*x[t])*(R+k*x[t]));
    J[t][i]=-2*b*(x[i]-x[j])/(R+k*x[t]);
    J[t][j]=2*b*(x[i]-x[j])/(R+k*x[t]);
    J[t][t]=1 + b*(x[i]-x[j])*(x[i]-x[j])*k/((R+k*x[t])*(R+k*x[t]));
}
void stamp2(MAT &J, VEC &x, int i, double V, int t) // special case for source and ground
{
    J[i][i]-=1/(R+k*x[t]);
    J[i][t]=(x[i]-V)*k/((R+k*x[t])*(R+k*x[t]));
    J[t][i]=-2*b*(x[i]-V)/(R+k*x[t]);
    J[t][t]=1 + b*(x[i]-V)*(x[i]-V)*k/((R+k*x[t])*(R+k*x[t]));
}
void J2(VEC &x, MAT &J, double V)   // jacobian for Q2
{
    stamp2(J, x, 0, V, 7);
    stamp(J, x, 0, 1, 8);
    stamp2(J, x, 2, V, 9);
    stamp(J, x, 0, 3, 10);
    stamp2(J, x, 1, gnd, 11);
    stamp(J, x, 2, 3, 12);
    stamp2(J, x, 3, gnd, 13);
    stamp(J, x, 2, 4, 14);
    stamp(J, x, 3, 5, 15);
    stamp2(J, x, 6, gnd, 16);
    stamp(J, x, 4, 5, 17);
    stamp(J, x, 5, 6, 18);
}
int main(void)  // main function
{
    fstream fout;
    fout.open("result.txt", fstream::out);  // write the data in txt file
    fout << "Voltage" << " " << "It" << " " << "I2" << " " << "I7" << " " << "I12" << endl;

    // The following code is for Q1
    /*
    double V;
    double err;
    for (V=0; V<=5; V+=0.1) {
        err=1+e;
        VEC x(7); VEC F(7); MAT J(7);
        F1(x, F, V);
        while (err>e) {
            VEC xt(7); VEC xtt(7); MAT T(7);
            J1(x, J, V);
            T=J;
            T=luFact(T);
            xt=fwdSubs(T, -F);
            xtt=bckSubs(T, xt);
            x += xtt;
            F1(x, F, V);
            err=norms(F, 2);
        }
        fout << V << " ";
        fout <<  (V-x[0])/(R+a*fabs(V-x[0])) + (V-x[2])/(R+a*fabs(V-x[2])) << " ";
        fout << (x[0]-x[1])/(R+a*fabs(x[0]-x[1]))  << " ";
        fout << (x[1])/(R+a*fabs(x[1])) << " ";
        fout << (x[5]-x[6])/(R+a*fabs(x[5]-x[6])) << endl;
    }
    */
    // Q2
    double V;
    double err;
    for (V=0; V<=5; V+=0.5) {
        err=1+e;
        VEC x(19); VEC F(19); MAT J(19);
        F2(x, F, V);
        while (err>e) {
            VEC xt(19); VEC xtt(19); MAT K(19);
            J2(x, J, V);
            K=J;
            K=luFact(K);
            xt=fwdSubs(K, F);
            xt=bckSubs(K, xt);
            x -= xt;
            F2(x, F, V);
            err=norms(F, 2);
            cout << err << endl;
        }
        fout << V << " ";
        fout <<  (V-x[0])/(R+a*fabs(V-x[0])) + (V-x[2])/(R+a*fabs(V-x[2])) << " ";
        fout << (x[0]-x[1])/(R+a*fabs(x[0]-x[1]))  << " ";
        fout << (x[1])/(R+a*fabs(x[1])) << " ";
        fout << (x[5]-x[6])/(R+a*fabs(x[5]-x[6])) << endl;
    }

    fout.close();
    return 0;
}
