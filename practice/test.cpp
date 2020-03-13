#include <iostream>
#include "MAT.h"
#include "VEC.h"
using namespace std;
int main()
{
    int *dim;
    int i, j, k;
    double temp;
    double s=0;
    cin >> *dim;
    MAT a(*dim);
    MAT b(*dim);
    for ( i=0; i< *dim; i++) {
        for ( j=0; j< *dim; j++){
            cin >> a[i][j];
        }
    }
    MAT T(a.dim());
    MAT TT(a.dim());
    MAT D(a.dim());
    for (int i=0; i<*dim; i++) {
        for (int j=0; j<*dim; j++){
            cout << a[i][j] << " ";
        }
        cout << '\n';
    }
    T[0]=a[0];
    for (k=1; k<*dim; k++) {
        for (i=0; i<k; i++) {
            temp=(a[k]*T[i])/(T[i]*T[i]);
            T[k]+=(T[i]*temp);
        }
        T[k]=a[k]-T[k];
    }
    TT=(T.tpose());
    D=T*TT;
    b=(D.tpose());
    for (i=0; i<*dim; i++) {
        for (j=0; j<*dim; j++){
            cout << b[i][j] << '\t';
        }
        cout << '\n';
    }
    for (i=0; i<*dim; i++) {
        for (j=0; (j<*dim) and (j!=i); j++) {
            s+=(b[i][j]*b[i][j]);
        }
    }
    cout << s << '\n';
    return 0;
}