/*  EE4070 Numerical Analysis
    HW06 Matrix Condition Number
    Name: Yu-Hsiu Huang     ID: 104061249   
    Date: 2018/04/12
*/
#include <iostream>
#include <sstream>
#include <cmath>
#include "MAT.h"
#include "VEC.h"
#define MAXITER 10000       // maximum iteration
#define TOL 0.000000001     // allowed tolerence
#define MU 0.000000001      // mu for shifted inverse power method
using namespace std;

int e=1;                    // 1 for e1, 2 for e2, 3 for e3, 4 for e4 
void Construct(MAT &A, int const num, double g, int source[2], int ground[2]);
void Symmetrize(MAT &A, int const num, double g, int source[2], int ground[2]);
int EVpwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);    // power method
int EViPwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);   // inverse power method
int EViPwrShift(MAT &A, VEC &q0, double &lambda, double mu, double tol, int maxiter);   // shifted inverse power method

int main(int argc, char *argv[]) 
{
    istringstream ss(argv[1]);
    int num;             // number of side-resistor
    ss >> num;
    int dim=(num+1)*(num+1);    // dimension of the matrix
    int ground[2];          // location of ground
    int source[2];          // location of source
    int step=0;
    double g=(double)num/2000;  // conductance
    double lambda_1;        // largest eigenvalue
    double lambda_n;        // smallest eigenvalue
    double kapa;            // condition number
    
    MAT A(dim);
    VEC q(dim);

    if (num<30) {
        source[0]=num/2;    // the location of node connected to the source
        source[1]=0;
        ground[0]=num;      // the location of grounded node
        ground[1]=num/2;
    }
    else {
        source[0]=0;        // the location of node connected to the source
        source[1]=num/2;
        ground[0]=num;      // the location of grounded node
        ground[1]=num/2;
    }

    q[0]=(double) 3/5;       // initial vector
    //q[1]=(double) 4/5;

    Construct(A, num, g, source, ground);
    Symmetrize(A, num, g, source, ground);
    step=EVpwr(A, q, lambda_1, TOL, MAXITER);
    step=EViPwr(A, q, lambda_n, TOL, MAXITER);
    //step=EViPwrShift(A, q, lambda_n, MU, TOL, MAXITER);
    kapa=lambda_1/lambda_n;
    cout << "largest eigenvalue: " << lambda_1 << '\n';
    cout << "smallest eigenvalue: " << lambda_n << '\n';
    cout << "Condition Number: " << kapa << '\n';
    return 0;
}
void Construct(MAT &A, int const num, double g, int source[2], int ground[2])    // construct and initialize the linear system
{
    // the verticle resistors
    for (int i=0; i < num; i++) {
        for (int j=0; j <= num; j++) {
            int n=i+j*(num+1);
            if (i==source[0] && j==source[1]) {
                A[n][n]=0.003;
                A[n+1][n]-=g;
                A[n+1][n+1]+=g;
            }
            else if (i+1==source[0] && j==source[1]) {
                A[n+1][n+1]=0.003;
                A[n][n+1]-=g;
                A[n][n]+=g;
            }
            else if (i==ground[0] && j==ground[1]) {
                A[n+1][n+1]=0.003;
                A[n][n+1]-=g;
                A[n][n]+=g;
            }
            else {
                A[n+1][n+1]+=g;
                A[n][n]+=g;
                A[n][n+1]-=g;
                A[n+1][n]-=g;
            }
        }
    }
    // the horizontal resistors
    for (int i=0; i < A.dim()-num-1; i++) {
        if (i==source[0]+source[1]*(num+1)) {
            A[i][i]=0.003;
            A[i+num+1][i]-=g;
            A[i+num+1][i+num+1]+=g; 
        }
        else if (i+num+1==source[0]+source[1]*(num+1)) {
            A[i+num+1][i+num+1]=0.003;
            A[i][i+num+1]-=g;
            A[i][i]+=g;
        }
        else if (i==ground[0]+ground[1]*(num+1)) {
            A[i][i]=0.003;
            A[i+num+1][i]-=g;
            A[i+num+1][i+num+1]+=g;
        }
        else if (i+num+1==ground[0]+ground[1]*(num+1)) {
            A[i+num+1][i+num+1]=0.003;
            A[i][i+num+1]-=g;
            A[i][i]+=g;
        }
        else {
            A[i][i]+=g;
            A[i+num+1][i+num+1]+=g;
            A[i][i+num+1]-=g;
            A[i+num+1][i]-=g;
        }
    }
}
void Symmetrize(MAT &A, int const num, double g, int source[2], int ground[2])   // symmetrize linear system
{
    int node_s=source[0]+source[1]*(num+1);     // node connected to source
    int node_g=ground[0]+ground[1]*(num+1);     // node connected to ground

    for (int i=0; i < A.dim(); i++) {
        if (i!=node_s && A[i][node_s]!=0) {
            A[i][node_s]=0;
        }
        else if (i!=node_g && A[i][node_g]!=0) {
            A[i][node_g]=0;
        }
        else;
    }
}
int EVpwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter)
{
    int step=0;
    double error=tol+1;
    double v=lambda;
    VEC z(q0.len());            
    VEC u(q0.len());
    VEC w(q0.len());
    VEC r(q0.len());
    VEC qk(q0.len());
    
    qk=q0;
    while(error>tol && step<=maxiter) {
        z=A*q0;
        q0=z/norms(z,2);
        lambda=q0*(A*q0);
        if (e==1) {         // error method 1
            error = fabs(lambda-v);
            v=lambda;
        }
        else if (e==2) {    // error method 1
            error = norms(q0-qk, 2);
            qk=q0;
        }
        else if (e==3) {    // error method 1
            r=A*q0-lambda*q0;
            error = norms(r, 2);
        }
        else {              // error method 1
            r=A*q0-lambda*q0;
            u=q0*A;
            w=u/norms(u,2);
            error = norms(r, 2)/(w*qk);
            qk=q0;
        }
        step++;
        //cout << step << " " << error << '\n';
    }
    return step;
}
int EViPwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter)
{
    int step=0;
    double error=tol+1;
    double v=lambda;
    MAT A0(A.dim());
    VEC z(q0.len());           
    VEC u(q0.len());
    VEC w(q0.len());
    VEC r(q0.len());
    VEC qk(q0.len());
    A0=A;
    A0=luFact(A0);
    while(error>tol && step<=maxiter) {
        qk=q0;
        z=fwdSubs(A0,qk);
        z=bckSubs(A0,z);
        q0=z/norms(z,2);
        lambda=q0*(A*q0);
        if (e==1) {         // error method 1
            error = fabs(lambda-v);
            v=lambda;
        }
        else if (e==2) {    // error method 1
            error = norms(q0-qk, 2);
        }
        else if (e==3) {    // error method 1
            r=A*q0-lambda*q0;
            error = norms(r, 2);
        }
        else {              // error method 1
            r=A*q0-lambda*q0;
            u=q0*A;
            w=u/norms(u,2);
            error = norms(r, 2)/(w*qk);
            qk=q0;
        }
        step++;
        //cout << step << " " << error << '\n';
    }
    return step;
}
int EViPwrShift(MAT &A, VEC &q0, double &lambda, double mu, double tol, int maxiter)
{
    int step=0;
    double error=tol+1;
    double v=lambda;
    MAT A0(A.dim());
    VEC z(q0.len());           
    VEC u(q0.len());
    VEC w(q0.len());
    VEC r(q0.len());
    VEC qk(q0.len());
    
    for (int i=0; i<A.dim(); i++) A0[i][i]=1;
    A0=A-mu*A0;
    A0=luFact(A0);
    while(error>tol && step<=maxiter) {
        qk=q0;
        z=fwdSubs(A0,qk);
        z=bckSubs(A0,z);
        q0=z/norms(z,2);
        lambda=q0*(A*q0);
        if (e==1) {         // error method 1
            error = fabs(lambda-v);
            v=lambda;
        }
        else if (e==2) {    // error method 1
            error = norms(q0-qk, 2);
        }
        else if (e==3) {    // error method 1
            r=A*q0-lambda*q0;
            error = norms(r, 2);
        }
        else {              // error method 1
            r=A*q0-lambda*q0;
            u=q0*A;
            w=u/norms(u,2);
            error = norms(r, 2)/(w*qk);
            qk=q0;
        }
        step++;
        //cout << step << " " << error << '\n';
    }
    return step;
}
