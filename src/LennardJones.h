#ifndef LENNARDJONES
#define LENNARDJONES

#include <armadillo>
#include "Structs.h"

using namespace arma;
using namespace QMC2;

/* Lennard-Jones potential helper functions */
void LJPotMatrix(GeneralParams &, mat&, const mat&);
double LJPotValue(GeneralParams &, const mat&);
double LJDiffValue(GeneralParams &, mat&, const mat&, const int&);
double LJPotTail(GeneralParams &);

/* Lennard-Jones virial helper functions */
void LJVirMatrix(GeneralParams &, mat&, const mat&);
double LJVirValueOverV(GeneralParams &, const mat&);
double LJUpdateVirValueOverV(GeneralParams &, mat&, const mat&, const int&);
double LJVirTail(GeneralParams &);

#endif // LENNARDJONES

