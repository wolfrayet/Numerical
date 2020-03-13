//  EE4070 Numerical Analysis
//  HW09 Spline Interpolations
//  Name: Yu-Hsiu Huang  ID: 104061249
//  Date: 2018/05/10
#include <iostream>
#include <math.h>
#include <fstream>
#include "VEC.h"
#include "MAT.h"
using namespace std;

void splineM(int N, VEC &X, VEC &Y, VEC &M);        // spline moment
double spline(double x, int N, VEC &X, VEC &Y, VEC &M);     // spline interpolate at x
int main(void) 
{
    fstream fin, fout;
    // ========== load data ============
    char xt, yt;            // only used for read in 
    double x_temp[21], y_temp[21];
    double x_real[301], y_real[301];
    double y = 0;           // interpolated value
    double e1 = 0;          // maximim error
    int n = 0;
    int i = 0;
    cin >> xt >> yt;
    cin >> x_temp[0] >> y_temp[0];
    while (i < 21 && x_temp[i] != 0) {
        i++;
        n++;                // calculate the dimension of the data
        cin >> x_temp[i] >> y_temp[i];
    }
    // ====== construct data vector =====
    VEC xdata(n);           // x from input file 
    VEC ydata(n);           // y from input file
    VEC M(n);               // boundary moment 
    for (int i=0; i<n; i++) {
        xdata[i]=x_temp[i];
        ydata[i]=y_temp[i];
    }
    // ========= read f301.dat ==========
    fin.open("f301.dat", fstream::in);
    fin >> xt >> yt;
    for (int i = 0; i < 301; i++) fin >> x_real[i] >> y_real[i];
    fin.close();
    // ============ txt file ============
    fout.open("result.txt",fstream::out);   // output in result.txt
    fout << "x"<< '\t' << "y" << '\n';
    // ========= interpolation ==========
    splineM(n, xdata, ydata, M);            // generate boundary moment
    for (double x = 475; x <= 775; x++) {   // find the interpolated y
        int i = (int) x;
        y = spline(x, n, xdata, ydata, M);
        fout << x << " " << y << '\n';      // write the result in result.txt
        if (fabs(y - y_real[i-475]) > e1) e1 = fabs(y - y_real[i-475]);
    }
    fout.close();
    cout << "No. of sample points: " << n << '\n';
    cout << "maximum error of interpolation: " << e1 << '\n';
    return 0;
}
void splineM(int N, VEC &X, VEC &Y, VEC &M) // generate boundary moment
{
    VEC d(N);   // need to solve the linear system: K*M = d to find M
    VEC t(N);   // temperal value for forward substution
    MAT K(N);
    for (int i=0; i<N; i++) {          
        for (int j=0; j<N; j++) {
            if (j == i) K[i][j]=2;
            else if ( i != 0 && i != N-1 && j == i-1) K[i][j]=(X[i]-X[i-1])/(X[i+1]-X[i-1]);
            else if ( i != 0 && i != N-1 && j == i+1) K[i][j]=(X[i+1]-X[i])/(X[i+1]-X[i-1]);
            else K[i][j] = 0;
        }
        if (i != 0 && i != N-1) {
            d[i]=6/(X[i+1]-X[i-1])*((Y[i+1]-Y[i])/(X[i+1]-X[i]) - (Y[i]-Y[i-1])/(X[i]-X[i-1]));
        }
        else d[i]=0;
    }
    // ======= LU decomposition ========
    K=luFact(K);
    t=fwdSubs(K, d);
    M=bckSubs(K, t);
}
double spline(double x, int N, VEC &X, VEC &Y, VEC &M)  // calculated interpolated y
{
    double y;
    int i=-1;       // x <= X[i] && x >= X[i-1]
    int h;          // lenth of the subinterval
    for (int index=0; index<N && i == -1; index++) {
        if (x == X[index]) {
            i=index;
            y = Y[index];
        }
        else if (x<X[index] && x>X[index-1]) {
            i=index;
            h=X[i]-X[i-1];
            y = M[i-1]/6/h*(X[i]-x)*(X[i]-x)*(X[i]-x) + M[i]/6/h*(x-X[i-1])*(x-X[i-1])*(x-X[i-1]);
            y = y + (Y[i-1] - M[i-1]/6*h*h)/h*(X[i]-x) + (Y[i] - M[i]/6*h*h)/h*(x-X[i-1]);
        }
    }
    return y;
}