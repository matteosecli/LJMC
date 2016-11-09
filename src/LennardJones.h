#ifndef LENNARDJONES
#define LENNARDJONES

#include <armadillo>
#include "Structs.h"

using namespace arma;
using namespace QMC2;

/* Lennard-Jones helper functions */
void LJPotMatrix(GeneralParams &, mat&, const mat&);
double LJPotValue(GeneralParams &, const mat&);
double LJDiffValue(GeneralParams &, mat&, const mat&, const int&);

#endif // LENNARDJONES

