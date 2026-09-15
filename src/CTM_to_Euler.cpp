#include "CTM_to_Euler.h"
#include <math.h>

bool CTM_to_Euler(const double C_b_n[9], double eul[3])
{
    // Column-major storage:
    // C(row,col) = C[row + 3*col]

    double s = -C_b_n[2];   // -C(3,1)

    if (s > 1.0) {
        s = 1.0;
    } else if (s < -1.0) {
        s = -1.0;
    }

    // roll
    eul[0] = atan2(C_b_n[5], C_b_n[8]);   // atan2(C(3,2), C(3,3))

    // pitch
    eul[1] = asin(s);                     // -asin(C(3,1))

    // yaw
    eul[2] = atan2(C_b_n[1], C_b_n[0]);   // atan2(C(2,1), C(1,1))

    return true;
}