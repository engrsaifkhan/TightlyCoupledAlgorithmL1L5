#include "Skew_Symmetric.h"

void Skew_Symmetric(const double a[3], double A[9])
{
    // Matrix:
    // [   0   -a3   a2
    //    a3    0   -a1
    //   -a2    a1   0 ]

    // Column-major storage:
    // [ A11 A21 A31 A12 A22 A32 A13 A23 A33 ]

    A[0] =  0.0;    // A11
    A[1] =  a[2];   // A21
    A[2] = -a[1];   // A31

    A[3] = -a[2];   // A12
    A[4] =  0.0;    // A22
    A[5] =  a[0];   // A32

    A[6] =  a[1];   // A13
    A[7] = -a[0];   // A23
    A[8] =  0.0;    // A33
}