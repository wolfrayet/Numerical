/*  MAT class functions
    Name: Yu-Hsiu Huang  ID: 104061249
    Date: 2018/04/12
*/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "MAT.h"
using namespace std;
MAT::MAT(int dim)               // uninit constructor
{
    n=dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++)
        va[i]=newVEC(n);
}
MAT::MAT(const MAT &m1)         // copy constructor
{
    VEC **vsrc=m1.va;           // to get around not indexing const MAT
    n=m1.n;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++) {
        va[i]=newVEC(n);
        (*va[i])=(*vsrc[i]);    // VEC assignment
    }
}
MAT::MAT(int dim, double *v)    // init constructor
{
    n=dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++) {
        va[i]=newVEC(n);
        for (int j=0; j<n; j++) 
            (*va[i])[j]=*(v++);     // array indexing + VEC indexing
    }
}
MAT::~MAT()                     // destructor
{
    for (int i=n-1; i>0; i--) {
        (*va[i]).~VEC();
    }
    free(va);
}
int MAT::dim()                       // return dimension of the matrix
{
    return n;
}
MAT MAT::tpose()                // matrix transpose
{
    MAT mnew(n);
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++) 
            mnew[i][j]=(*va[j])[i];
    return mnew;
}
MAT &MAT::operator-()           // unary operator, negative value
{
    for (int i=0; i<n; i++) 
        for (int j=0; j<n; j++) 
            (*va[i])[j]=-(*va[i])[j];
    return *this;
}
MAT &MAT::operator=(MAT m1)      // assignment
{
    for (int i=0; i<n; i++)     // VEC assignment
        (*va[i])=m1[i];
    return *this;
}
MAT &MAT::operator+=(MAT &m1)   // m += m1
{
    for (int i=0; i<n; i++) 
        (*va[i])+=m1[i];        // VEC += operator
    return *this;
}
MAT &MAT::operator-=(MAT &m1)   // m -= m1
{
    for (int i=0; i<n; i++) 
        (*va[i])-=m1[i];        // VEC -= operator
    return *this;
}
MAT &MAT::operator*=(double db1)    // m *= db1
{
    for (int i=0; i<n; i++)
        (*va[i])*=db1;          // VEC *= operator
    return *this;
}
MAT &MAT::operator/=(double db1)    // m /= db1
{
    for (int i=0; i<n; i++)
        (*va[i])/=db1;          // VEC /= operator
    return *this;
}
MAT MAT::operator+(MAT m1)      // matrix addition
{
    MAT z(n);
    for (int i=0; i<n; i++)
        z[i]=(*va[i])+m1[i];    // VEC addition and assignment
    return z;
}
MAT MAT::operator-(MAT m1)      // matrix subtraction
{
    MAT z(n);
    for (int i=0; i<n; i++)
        z[i]=(*va[i])-m1[i];    // VEC subtraction and assignment
    return z;
}
MAT MAT::operator*(MAT m1)      // matrix multplication
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=0;
            for (int k=0; k<n; k++)
                z[i][j]+=((*va[i])[k]*m1[k][j]);
        }
    }
    return z;
}
VEC &MAT::operator[](int m)     // m'th row of the matrix
{                
    return (*va[m]);
}
VEC MAT::operator*(VEC v1)      // m * v1
{
    VEC s(n);
    for (int i=0; i<n; i++) 
        s[i]=(*va[i])*v1;       // VEC inner product
    return s;
}
MAT MAT::operator*(double db1)  // Matrix multiplied by scalar
{
    MAT z(n);
    for (int i=0; i<n; i++) 
        z[i]=((*va[i])*db1);    // VEC multiplied by scalar
    return z;
}
MAT MAT::operator/(double db1)  // Matrix divided by scalar
{
    MAT z(n);
    for (int i=0; i<n; i++)
        z[i]=((*va[i])/db1);    // VEC divided by scalar
    return z;
}
MAT operator*(double db1, const MAT &m1)    // db1 x m
{
    MAT z(m1.n);
    for (int i=0; i<m1.n; i++)
        z[i]=(*m1.va[i]*db1);   // VEC multiplied by scalar
    return z;
}
VEC operator*(VEC &v1, MAT &m1) // vT x m
{
    VEC v2(m1.n);
    for (int i=0; i<m1.n; i++) {
        v2[i]=0;
        for (int j=0; j<m1.n; j++)
            v2[i] += v1[j]*m1[j][i];
    }
    return v2;
}
MAT &luFact(MAT &m1)            // LU decomposition 
{
    int i, j, k;
    for (i=0; i<m1.n; i++) {
        for (j=i+1; j<m1.n; j++) {
            m1[j][i]/=m1[i][i]; // form l[j][i]
        }
        for (j=i+1; j<m1.n; j++) {          // update lower matrix
            for (k=i+1; k<m1.n; k++) {
                m1[j][k]-=m1[j][i]*m1[i][k];
            }
        }
    }
    return m1;                  // the return matrix = L + U - I
}
VEC fwdSubs(MAT &m1, VEC b)     // forward substitution
{
    VEC y(m1.n);
    int i, j;
    for (i=0; i<m1.n; i++) {
        y[i]=b[i];              // initialize y to b
    }
    for (i=0; i<m1.n; i++) {    // solve L * y = b
        for (j=i+1; j<m1.n; j++) {
            y[j]-=m1[j][i]*y[i];
        }
    }
    return y;                   // return root y
}
VEC bckSubs(MAT &m1, VEC y)     // backword substitution
{
    VEC x(m1.n);
    int i, j;
    for (i=0; i<m1.n; i++) {
        x[i]=y[i];              // initialize x to y
    }
    for (i=m1.n-1; i>=0; i--) { // solve U * x = y
        x[i]/=m1[i][i];
        for (j=i-1; j>=0; j--) {
            x[j]-=m1[j][i]*x[i];
        }
    }
    return x;                   // return root x
}
int jacobi(MAT &A, VEC b, VEC &x, int maxIter, double tol)      // Jacobi method
{
    VEC xk(x.len());            // x(k) term
    double error=1;
    double t1;                  // temporary value
    int p=0;                    // p=0 for infinite norm
    int step=0;
    while (error > tol) {
        xk=x;
        for (int i=0; i<A.dim(); i++) {
            t1=0;
            if (i!=10 && i!=230) {
                for (int j=0; j<A.dim(); j++) {
                    if (j!=i) t1+=A[i][j]*xk[j];
                }
                x[i]=(b[i]-t1)/A[i][i];
            }
        }
        error=norms(x-xk, p);
        step+=1;
    }
    if (step <= maxIter) return step;
    else return maxIter;
}
int gaussSeidel(MAT &A, VEC b, VEC &x, int maxIter, double tol)     // Gauss-Seidel method
{
    VEC xk(x.len());            // x(k) term
    double error=1;
    double t1=0;                // temporary values
    double t2=0;
    int p=0;                    // p=0 for infinite norm
    int step=0;
    while (error > tol) {
        xk=x;
        for (int i=0; i<A.dim(); i++) {     // forward substitution
            for (int j=0; j<=i-1 ; j++) {
                t1+=A[i][j]*x[j];
            }
            for (int j=i+1; j<A.dim(); j++) {
                t2+=A[i][j]*xk[j];
            }
            x[i]=(b[i]-t1-t2)/A[i][i];
            t1=0;
            t2=0;
        }        
        error=norms(x-xk, p);
        step+=1;
    }
    if (step <= maxIter) return step;
    else return maxIter;
}
int sgs(MAT &A, VEC b, VEC &x, int maxIter, double tol)     // Symmetric Gauss-Seidel method
{
    VEC xk(x.len());            // x(k) term
    VEC xkk(x.len());           // x(k+1/2) term
    double error=1;
    double t1=0;                // temporary values
    double t2=0;
    int p=0;                    // p=0 for infinite norm
    int step=0;
    while (error > tol) {
        xk=x;
        for (int i=0; i<A.dim(); i++) {     // forward substutution
            t1=0;
            t2=0;
            for (int j=0; j<=i-1 ; j++) {
                t1+=A[i][j]*xkk[j];
            }
            for (int j=i+1; j<A.dim(); j++) {
                t2+=A[i][j]*xk[j];
            }
            xkk[i]=(b[i]-t1-t2)/A[i][i];
        }
        for (int i=A.dim()-1; i>=0; i--) {  // backward substitution
            t1=0;
            t2=0;
            for (int j=0; j<=i-1 ; j++) {
                t1+=A[i][j]*xkk[j];
            }
            for (int j=i+1; j<A.dim(); j++) {
                t2+=A[i][j]*x[j];
            }
            x[i]=(b[i]-t1-t2)/A[i][i];
        }        
        error=norms(x-xk, p);   // calculate norm of error
        step+=1;
    }
    if (step <= maxIter) return step;
    else return maxIter;
}
int cg(MAT &A, VEC b, VEC &x, int maxIter, double tol)  // Conjugate Gradient Method
{
    VEC xk(x.len());
    VEC rk(x.len());
    VEC pk(x.len());
    double alpha=0;
    double beta=0;
    double error=1;
    int step=0;
    rk=b;
    pk=rk;
    while (error > tol) {                 // iteration
        xk=x;
        alpha=(pk*rk)/(pk*(A*pk));
        x=xk+alpha*pk;
        rk-=alpha*(A*pk);
        beta=(pk*(A*rk))/(pk*(A*pk));
        pk=rk-beta*pk;
        error=sqrt((rk*rk)/x.len());
        step+=1;
    }
    if (step <= maxIter) return step;
    else return maxIter;
}
int cg2(MAT &A, VEC b, VEC &x, int maxIter, double tol)  // Conjugate Gradient Method, 2nd form
{
    VEC xk(x.len());
    VEC rk(x.len());
    VEC rkk(x.len());
    VEC pk(x.len());
    double alpha=0;
    double beta=0;
    double error=1;
    int step=0;
    rk=b;               
    pk=b;
    while (error > tol) {               // iteration
        xk=x;
        alpha=(rk*rk)/(pk*(A*pk));
        x=xk+alpha*pk;
        rkk=rk-alpha*(A*pk);
        beta=(rkk*rkk)/(rk*rk);
        pk=rkk+beta*pk;
        rk=rkk;
        error=sqrt((rk*rk)/x.len());
        step+=1;
    }
    if (step <= maxIter) return step;
    else return maxIter;
}
void QR(MAT &A, MAT &Q, MAT &R)         // QR factorization
{
    A = A.tpose();
    R[0][0]=sqrt(A[0]*A[0]);
    Q[0]=A[0]/R[0][0];
    for (int j=1; j<A.dim(); j++) {
        Q[j]=A[j];
        for (int i=0; i<j; i++) {
            R[i][j]=Q[i]*Q[j];
            Q[j]=Q[j]-R[i][j]*Q[i];
        }
        R[j][j]=sqrt(Q[j]*Q[j]);
        Q[j]=Q[j]/R[j][j];
    }
    Q = Q.tpose();
}
int EVqr(MAT &A, double tol, int maxiter)   // QR iteration
{
    int step=0;
    double error=0.1;
    MAT Q(A.dim());
    MAT R(A.dim());
    while (error > tol && step < maxiter) {
        error=0;                            // initialize
        // ========iteration==========
        QR(A, Q, R);
        A = R*Q;
        // ==========error============
        for (int i=1; i<A.dim(); i++) {
            if (fabs(A[i][i-1]) > error) error=fabs(A[i][i-1]);
        }
        // ==========steps============
        step++;
    }
    if (step < maxiter) return step;
    else return maxiter;
}
int EVqrshifted(MAT &A, double mu, double tol, int maxiter) // QR iteration with shifted
{
    int step=0;
    double error=tol+1;
    MAT I(A.dim());
    MAT Q(A.dim());
    MAT R(A.dim());

    for (int i=0; i<I.dim(); i++) I[i][i]=mu;   // initialize unit matrix
    while (error > tol && step < maxiter) {
        error=0;
        // ========iteration==========
        A -= I;
        QR(A, Q, R);
        A = R*Q + I;
        // ==========error============
        for (int i=1; i<A.dim(); i++) {
            if (fabs(A[i][i-1]) > error) error=fabs(A[i][i-1]);
        }
        // ==========steps============
        step++;
    }
    if (step < maxiter) return step;
    else return maxiter;
}