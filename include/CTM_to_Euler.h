#ifndef CTM_TO_EULER_H
#define CTM_TO_EULER_H

#include <math.h>
#include <stdbool.h>

// Keeps the same interface as your current project
bool CTM_to_Euler(const double C[9], double eul[3]);

#endif  // EULER_TO_CTM_H