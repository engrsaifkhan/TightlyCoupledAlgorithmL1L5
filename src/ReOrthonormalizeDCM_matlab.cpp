//
// Optimized ReOrthonormalizeDCM_matlab.cpp
//
// Optimized SVD-based DCM re-orthonormalization for embedded INS/TCA use.
// Keeps SVD approach but removes unnecessary MATLAB Coder runtime overhead.
//

#include "ReOrthonormalizeDCM_matlab.h"
#include "svd.h"
#include <cmath>

namespace
{
    // Optional determinant calculation for 3x3 DCM
    double det3x3(const double C[9])
    {
        return C[0] * (C[4] * C[8] - C[5] * C[7])
             - C[3] * (C[1] * C[8] - C[2] * C[7])
             + C[6] * (C[1] * C[5] - C[2] * C[4]);
    }
}


//
// Arguments:
//   C_b_e : column-major 3x3 DCM
//
// Operation:
//   C = U*S*V^T
//   C_ortho = U*V^T
//
//
void ReOrthonormalizeDCM_matlab(double C_b_e[9])
{
    // Static workspace avoids repeated stack allocation
    static double U[9];
    static double V[9];
    static double S[3];

    bool valid = true;

    for (int i = 0; i < 9; i++)
    {
        if (!std::isfinite(C_b_e[i]))
        {
            valid = false;
            break;
        }
    }

    if (!valid)
    {
        return;
    }


    //
    // Singular Value Decomposition
    //
    // C = U*S*V'
    //
    coder::internal::svd(C_b_e, U, S, V);


    //
    // Polar correction:
    //
    // C_ortho = U * V^T
    //
    double C_new[9];

    for (int col = 0; col < 3; col++)
    {
        for (int row = 0; row < 3; row++)
        {
            C_new[row + 3 * col] =
                U[row]
                * V[col]
                +
                U[row + 3]
                * V[col + 3]
                +
                U[row + 6]
                * V[col + 6];
        }
    }


    //
    // Ensure proper rotation matrix:
    //
    // det(C)=+1
    //
    if (det3x3(C_new) < 0.0)
    {
        // Flip third column of U
        for (int row = 0; row < 3; row++)
        {
            U[row + 6] = -U[row + 6];
        }


        for (int col = 0; col < 3; col++)
        {
            for (int row = 0; row < 3; row++)
            {
                C_new[row + 3 * col] =
                    U[row]
                    * V[col]
                    +
                    U[row + 3]
                    * V[col + 3]
                    +
                    U[row + 6]
                    * V[col + 6];
            }
        }
    }


    //
    // Copy result back
    //
    for (int i = 0; i < 9; i++)
    {
        C_b_e[i] = C_new[i];
    }
}
