/*  EE4070 Numerical Analysis
    HW3. Resistor Networks
    Name: Yu-Hsiu Huang     ID: 104061249
    Date: 2018/03/23
*/
#include <iostream>
#include "MAT.h"
#include "VEC.h"
using namespace std;

double Resistor(int const Num);     // assign the conductance
void Construct(MAT &A, VEC &b, int const Num, double g);    // construct the linear system
void EquiRes(VEC &x, int const Num, double g);      // find the equivalent resistance
void GetAnswer(VEC &x, int const Num);      // output the desired answer

int main(void)
{
    int Num;                // number of resistors each side
    int dim;                // dimension of linear system
    double g;               // conductance of each resistor
    Num=50;                 // change the value of Num to switch the simulation mode
    dim=(Num+1)*(Num+1);    // dimension of the linear system
    MAT A(dim);             // matrix for node calculation
    VEC b(dim);             // right-hand side
    VEC y(dim);
    VEC x(dim);             // node voltage

    g=Resistor(Num);        // assign the conductance value
    Construct(A, b, Num, g);    // construct the linear system
    A=luFact(A);            // LU decomposition
    y=fwdSubs(A, b);        // forward substitution
    x=bckSubs(A, y)/g;      // backward substitution
    EquiRes(x, Num, g);     // find equivalent resistor
    GetAnswer(x, Num);      // print the answers
    return 0;
}
double Resistor(int const Num)      // assign the conductance value 
{
    if (Num==2) return 0.001;
    else if (Num==4) return 0.002;
    else if (Num==10) return 0.005;
    else if (Num==20) return 0.01;
    else if (Num==40) return 0.02;
    else return 0.025;
}
void Construct(MAT &A, VEC &b, int const Num, double g)     // construct the linear system
{
    // the right-hand side of linear system
    for (int i=0; i<b.len(); i++) { 
        if (i!=(Num/2)) b[i]=0;
        else b[i]=1;
    }
    // construct matrix through verticle resistors
    for (int i=0; i<Num; i++) { 
        for (int j=0; j<=Num; j++) {
            int n=i+j*(Num+1);
            if (i==Num/2 && j==0) {
                A[n][n]=1/g;
                A[n+1][n]-=1;
                A[n+1][n+1]+=1;
            }
            else if (i+1==Num/2 && j==0) {
                A[n+1][n+1]=1/g;
                A[n][n+1]-=1;
                A[n][n]+=1;
            }
            else if (i==Num-1 && j==Num/2) {
                A[n+1][n+1]=1/g;
                A[n][n+1]-=1;
                A[n][n]+=1;
            }
            else {
                A[n+1][n+1]+=1;
                A[n][n]+=1;
                A[n][n+1]-=1;
                A[n+1][n]-=1;
            }
        }
    }
    // construct matrix through horizontal resistors
    for (int i=0; i<b.len()-Num-1; i++) {
        if (i==Num/2) {
            A[i][i]=1/g;
            A[i+Num+1][i]-=1;
            A[i+Num+1][i+Num+1]+=1;
        }
        else if (i==(Num/2+1)*(Num+1)-1) {
            A[i][i]=1/g;
            A[i+Num+1][i]-=1;
            A[i+Num+1][i+Num+1]+=1;
        }
        else if (i+Num+1==(Num/2+1)*(Num+1)-1) {
            A[i+Num+1][i+Num+1]=1/g;
            A[i][i+Num+1]-=1;
            A[i][i]+=1;
        }
        else {
            A[i][i]+=1;
            A[i+Num+1][i+Num+1]+=1;
            A[i][i+Num+1]-=1;
            A[i+Num+1][i]-=1;
        }
    }
}
void EquiRes(VEC &x, int const Num, double g)   // find the equivalent resistor
{
    double current;     // net current through the circuit
    current=g*(3*x[Num/2]-x[Num/2+1]-x[Num/2-1]-x[Num/2+Num+1]);    // apply Kirchhoff's current law
    cout << "Equivalent Resistance: " << 1/current << '\n'; 
}
void GetAnswer(VEC &x, int const Num)   // output the desired answer
{
    cout << "Voltage value for the south-west corner: " << x[Num] << '\n';
    cout << "Voltage value for the north-east corner: " << x[x.len()-Num-1] << '\n';
    cout << "Voltage value for the center of east side: " << x[x.len()-Num/2-1] << '\n';
}
