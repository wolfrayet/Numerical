// matrix class
#ifndef MAT_H
#define MAT_H
#include "VEC.h"

class MAT 
{
    private:
        int n;                  // define nxn matrix
        VEC **va;               // array of n pointers to vectors
    public:
        MAT(int dim);           // uninit constructor
        MAT(const MAT &m1);     // copy constructor
        MAT(int dim, double *v);    // init constructor
        ~MAT();                 // destructor
        int dim();              // return dimension of the matrix
        MAT tpose();            // transpose
        MAT &operator-();       // unary operator, negative value
        MAT &operator=(MAT m1);     // assignment
        MAT &operator+=(MAT &m1);   // m += m1
        MAT &operator-=(MAT &m1);   // m -= m1
        MAT &operator*=(double db1);    // m *= db1
        MAT &operator/=(double db1);    // m /= db1
        MAT operator+(MAT m1);  // m + m1
        MAT operator-(MAT m1);  // m - m1
        MAT operator*(MAT m1);  // m * m1
        VEC &operator[](int m);    // m'th row
        VEC operator*(VEC v1);  // m X v1
        MAT operator*(double db1);  // m * db1
        MAT operator/(double db1);  // m / db1
        friend MAT operator*(double db1, const MAT &m1);  // db1 x m
        friend VEC operator*(VEC &v1, MAT &m1);     // vT x m
};
MAT operator*(double db1, const MAT &m1);   // db1 x m
VEC operator*(VEC &v1, MAT &m1);    // vT x m
#endif
