//  EE4070 Numerical Analysis
//  HW08 Polynomial Interpolations
//  Name: Yu-Hsiu Huang ID: 104061249
//  Date: 2018/05/03
#include <iostream>
#include <math.h>
#include <fstream>
#include "VEC.h"
#include "MAT.h"
using namespace std;

double Lagrange(double x, VEC &XDATA, VEC &YDATA);
int main(void)
{
    fstream fin, fout;
    // ======= load the data ==========
    char xt, yt;
    double x_temp[21], y_temp[21];
    double x_real[301], y_real[301];
    double y = 0;
    double e1 = 0;
    double e2 = 0;
    int n = 0;
    int i = 0;
    cin >> xt >> yt;
    cin >> x_temp[0] >> y_temp[0];
    while (i < 21 && x_temp[i] != 0) {
        i++;
        n++;
        cin >> x_temp[i] >> y_temp[i];
    }
    VEC xdata(n);
    VEC ydata(n);
    for (int i = 0; i < n && x_temp[i] != 0; i++) {
        xdata[i] = x_temp[i];
        ydata[i] = y_temp[i];
    }
    // ======== read f301.dat ===========
    fin.open("f301.dat", fstream::in);
    fin >> xt >> yt;
    for (int i = 0; i < 301; i++) {
        fin >> x_real[i] >> y_real[i];
    }
    fin.close();
    // ======== output ==================
    cout << "No. of sample points: " << n << '\n';
    fout.open("result.txt",fstream::out);
    fout << "x"<< '\t' << "y" << '\n';
    // ======= lagrange interpolation ===
    for (double x = 475; x <= 775; x++) {
        int i = (int)x;
        y = Lagrange(x, xdata, ydata);
        fout << x << " " << y << '\n';
        if (fabs(y - y_real[i-475]) > e1) { e1 = fabs(y - y_real[i-475]); n = i;}
        if (i >= 550 && i <= 700 && fabs(y - y_real[i-475]) > e2) e2 = fabs(y - y_real[i-475]);
    }
    fout.close();
    cout << "maximum error of interpolation: " << e1  << " " << n << '\n';
    cout << "maximum error in range [550,700]: " << e2 << '\n';
    return 0;
}
double Lagrange(double x, VEC &XDATA, VEC &YDATA)
{
    VEC y(XDATA.len());
    for (int i = 0; i < XDATA.len(); i++) y[i] = YDATA[i];
    for (int k = 1; k < XDATA.len(); k++) {
        for (int j = 0; j < XDATA.len() - k; j++) {
            y[j] = ((x - XDATA[j])*y[j+1] - (x - XDATA[k+j])*y[j])/(XDATA[j+k] - XDATA[j]);
        }
    }
    return y[0];
}