// vector class
#ifndef VEC_H
#define VEC_H

class VEC 
{
    private:
        int dim;                    // vector dimension 
        double *val;                // array to store vector
    public:                     
        VEC(int n);                 // uninit constructor, val set to 0
        VEC(const VEC &v1);         // copy constructor
        VEC(int n, double *v);      // init constructor
        ~VEC();                     // destructor
        int len();                  // dimension of the vector
        VEC &operator-();           // unitary operator, negative value
        VEC &operator=(const VEC v1);   // assignment
        VEC &operator+=(const VEC v1);  // V += v1;
        VEC &operator-=(const VEC v1);  // V -= v1;
        VEC &operator*=(double db1);    // V *= db1;
        VEC &operator/=(double db1);    // V /= db1;
        VEC operator+(const VEC v1);    // V + v1
        VEC operator-(const VEC v1);    // V - v1
        double operator*(VEC v1);   // inner product
        VEC operator*(double db1);  // V * db1
        VEC operator/(double db1);  // V / db1
        double &operator[](int n);  // indexing
        friend VEC operator*(double db1, const VEC v1);     // db1 x V
        friend VEC *newVEC(int n);  // create dynamic VEC
};
VEC operator*(double db1, const VEC v1);
VEC *newVEC(int n);                 // create dynamic VEC
#endif
