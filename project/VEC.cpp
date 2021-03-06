// VEC class functions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VEC.h"

VEC::VEC(int n)             // uninit constructor                           
{
    dim=n;
    val=(double *)calloc(n, sizeof(double));
}
VEC::VEC(const VEC &v1)     // copy constructor
{
    dim=v1.dim;
    val=(double *)calloc(dim, sizeof(double));
    for (int i=0; i<dim; i++) {
        val[i]=v1.val[i];
    }
}
VEC::VEC(int n, double *v)  // init constructor
{
    dim=n;
    val=(double *)calloc(n, sizeof(double));
    for (int i=0; i<n; i++) val[i]=v[i];
}
VEC::~VEC()                 // destructor
{                          
    free(val);
}
int VEC::len()              // return dimension of the vector
{
    return dim;
}
VEC &VEC::operator-()       // unary operator - : negative value
{
    for (int i=0; i<dim; i++) {
        val[i]=-val[i];
    }
    return *this;
}
VEC &VEC::operator=(const VEC v1)   // assignment
{
    dim=v1.dim;
    for (int i=0; i<dim; i++) {
        val[i]=v1.val[i];
    }
    return *this;
}
VEC &VEC::operator+=(const VEC v1)  // V += v1
{
    for (int i=0; i<dim; i++) {
        val[i]+=v1.val[i];
    }
    return *this;
}
VEC &VEC::operator-=(const VEC v1)  // V -= v1
{
    for (int i=0; i<dim; i++) {
        val[i]-=v1.val[i];
    }
    return *this;
}
VEC &VEC::operator*=(double db1)    // V *= db1
{
    for (int i=0; i<dim; i++) {
        val[i]*=db1;
    }
    return *this;
}
VEC &VEC::operator/=(double db1)    // V /= db1
{
    for (int i=0; i<dim; i++) {
        val[i]/=db1;
    }
    return *this;
}
VEC VEC::operator+(const VEC v1)    // V + v1
{
    VEC s(*this);
    for (int i=0; i<dim; i++) {
        s.val[i]+=v1.val[i];
    }
    return s;
}
VEC VEC::operator-(const VEC v1)    // V - v1
{
    VEC s(*this);
    for (int i=0; i<dim; i++) {
        s.val[i]-=v1.val[i];
    }
    return s;
}
double VEC::operator*(VEC v1)       // inner product
{
    double p=0;
    for (int i=0; i<dim; i++) {
        p+=(val[i]*v1.val[i]);
    }
    return p;
}
VEC VEC::operator*(double db1)      // V * db1
{
    VEC s(*this);
    for (int i=0; i<dim; i++) {
        s.val[i]*=db1;
    }
    return s;
}
VEC VEC::operator/(double db1)      // V * db1
{
    VEC s(*this);
    for (int i=0; i<dim; i++) {
        s.val[i]/=db1;
    }
    return s;
}
double &VEC::operator[](int n)      // indexing
{
    if (n<0) n=0;
    else if (n>=dim) n=dim-1;
    return val[n];
}
VEC operator*(double db1, const VEC v1)     // db1 x v1
{
    VEC s(v1.dim, v1.val);
    for (int i=0; i<v1.dim; i++) 
        s.val[i]=(v1.val[i]*db1);   // V * db1
    return s;
}
VEC *newVEC(int n)                  // allocate dynamic VEC
{
    VEC *vptr;
    vptr=(VEC *)malloc(sizeof(VEC));
    vptr->dim=n;
    vptr->val=(double *)calloc(n, sizeof(double));
    return vptr;
}
double norms(VEC x, int p)          // calculate p-norm of a vector
{
    double norm=0;
    if (p==1) {
        for (int i=0; i<x.len(); i++) {
            norm+=abs(x[i]);
        }
    }
    else if (p==2) {
        for (int i=0; i<x.len(); i++) {
            norm+=x[i]*x[i];
        }
        norm=sqrt(norm);
    }
    else {                          // p=0 stands for infinity norm
        norm=x[0];
        for (int i=1; i<x.len(); i++) {
            if (x[i]>norm) norm=x[i];
        }
    }
    return norm;
}
