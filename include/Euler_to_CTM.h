#ifndef EULER_TO_CTM_H
#define EULER_TO_CTM_H

#include <math.h>
#include <stdbool.h>

// Keeps the same interface as your current project
bool Euler_to_CTM(const double eul[3], double C[9]);

#endif  // EULER_TO_CTM_H