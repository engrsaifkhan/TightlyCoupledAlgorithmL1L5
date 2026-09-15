#ifndef SKEW_SYMMETRIC_H
#define SKEW_SYMMETRIC_H

// Input:
//   a[3] = [a1, a2, a3]
//
// Output:
//   A[9] = 3x3 skew-symmetric matrix stored in column-major order
//          [A11 A21 A31 A12 A22 A32 A13 A23 A33]
void Skew_Symmetric(const double a[3], double A[9]);

#endif