//  EE4070 Numerical Analysis
//  Project. Parametric Spline Interpolation
//  Yu-Hsiu Huang 104061249
//  15 June, 2018
#include <iostream>
#include <math.h>
#include <fstream>
#include "VEC.h"
#include "MAT.h"
#define e 0.000001
#define MAX 1000
using namespace std;

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
double find(double x, int N, VEC &T, VEC &X, VEC &M) // find the corresponding t for x
{
    double t=1;
    int k=0;
    int n=3;
    double a[4], b[4], c[3], z[3];
    int find=0;
    for (int i=0; i<N && find==0; i++) {
        if (x<X[i] && x>X[i-1]) {
            find=1;
            k=i;
        }
        else if (x==X[i]) {
            find=2;
            t=T[i];
        }
    }
    double h=T[k]-T[k-1];
    // coeffiecient of 3rd polynomial
    a[0]=X[k-1]-x;
    a[1]=(X[k]-X[k-1])/h-h*(M[k]+2*M[k-1])/6;
    a[2]=M[k-1]/2;
    a[3]=(M[k]-M[k-1])/6/h;
    while (n>=1 && find==1) {               // solve the 3rd order polynomial to find t
        double err=1+e;
        int step=0;
        while (err>=e && step<MAX) {
            b[n]=a[n]; c[n-1]=b[n];
            for (int j=n-1; j>=0; j--) b[j]=a[j]+t*b[j+1];
            for (int j=n-2; j>=0; j--) c[j]=b[j+1]+t*c[j+1];
            t=t-b[0]/c[0];
            err=fabs(b[0]); step++;
        }
        z[n-1]=t;
        for (int j=0; j<n; j++) a[j]=b[j+1];
        n--;
    }
    for (int i=0; i<3 && find==1; i++) {
        if (z[i] < 1 && z[i]>0) {
            t=z[i]+T[k-1];
        }
    }
    return t;
}
int main(void) 
{
    fstream fin, fout;
    // ========== load data ============
    double x_real[301], y_real[301];
    double r;
    double y = 0;           // interpolated value
    double e1 = 0;          // maximim error
    int n = 0;
    cin >> n;
    // ====== construct data vector =====
    VEC xdata(n);           // x from input file 
    VEC ydata(n);           // y from input file
    VEC t(n);               // parameter t
    VEC Mx(n);               // boundary moment 
    VEC My(n);
    for (int i=0; i<n; i++) {
        cin >> xdata[i] >> ydata[i];
        t[i]=i;
    }
    // ========= read f301.dat ==========
    fin.open("f301.dat", fstream::in);
    fin >> r;
    for (int i = 0; i < 301; i++) fin >> x_real[i] >> y_real[i];
    fin.close();
    // ============ txt file ============
    fout.open("result.txt",fstream::out);   // output in result.txt
    fout << "t x y error" << '\n';
    // ========= interpolation ==========
    splineM(n, t, xdata, Mx);
    splineM(n, t, ydata, My);            // generate boundary moment
    for (double x = 475; x <= 775; x++) {   // find the interpolated y
        int i = (int) x;
        double par=find(x, n, t, xdata, Mx);
        y = spline(par, n, t, ydata, My);
        fout << par << " " << x << " " << y << " " << fabs(y - y_real[i-475]) << '\n';      // write the result in result.txt
        if (fabs(y - y_real[i-475]) > e1) e1 = fabs(y - y_real[i-475]);
    }
    fout.close();
    cout << "No. of sample points: " << n << '\n';
    cout << "maximum error of interpolation: " << e1 << '\n';
    return 0;
}
