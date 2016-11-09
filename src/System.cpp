#include "System.h"

#include <math.h>
#include <sstream>
#include <iomanip>

#include "Potential.h"


using namespace QMC2;

System::System(int number_particles, int dimension) {
    this->n_p = number_particles;
    this->dim = dimension;
}

System::System() {};

void System::add_potential(Potential* pot) {
    potentials.push_back(pot);
}

double System::get_potential_energy(const mat& r) {
    double potE = 0;
    double potE_i;

    for (Potential *pot : potentials)
    {
        potE_i = pot->get_pot_E(r);
        potE += potE_i;
    }

    return potE;
}
