#include <armadillo>
#include "Potential.h"
#include "Structs.h"

using namespace arma;
using namespace QMC2;

// Just the constructor
Potential::Potential(int number_particles, int dimension, std::string name) {
    this->n_p = number_particles;
    this->dim = dimension;
    this->name = name;
};

Potential::Potential () {};


// Harmonic oscillator definitions
Harmonic_osc::Harmonic_osc( GeneralParams & gP ) :
    Potential ( gP.number_particles, gP.dimension, "Harmonic Oscillator") {
        this->w = gP.frequency;
};

double Harmonic_osc::get_pot_E(const mat& r) const {
    int i = 0, j = 0;
    double e_potential = 0, r_single_particle = 0;
    for (i = 0; i < n_p; i++) {
      r_single_particle = 0;
      for (j = 0; j < dim; j++) {
        r_single_particle += pow(r(i,j),2);
      }
      e_potential += 0.5*pow(w,2)*r_single_particle;
    }
    return e_potential;
}


// Coulomb definitions
Coulomb::Coulomb( GeneralParams & gP ) :
    Potential ( gP.number_particles, gP.dimension, "Nuclear Coulomb Attraction") {
        this->Z = gP.nuclear_charge;
};

double Coulomb::get_pot_E(const mat& r) const {
    int i = 0, j = 0;
    double e_potential = 0, r_single_particle = 0;
    for (i = 0; i < n_p; i++) {
      r_single_particle = 0;
      for (j = 0; j < dim; j++) {
        r_single_particle += pow(r(i,j),2);
      }
      e_potential -= Z/r_single_particle;
    }
    return e_potential;
}


// eRepulsion definitions
eRepulsion::eRepulsion(GeneralParams & gP) :
    Potential ( gP.number_particles, gP.dimension, "Electron Repulsion") {
};

double eRepulsion::get_pot_E(const mat& r) const {
    int i = 0, j = 0, k = 0;
    double e_potential = 0, r_12 = 0;
    for (i = 0; i < n_p-1; i++) {
        for (j = i+1; j < n_p; j++) {
            r_12 = 0;
            for (k = 0; k < dim; k++) {
                r_12 += pow((r(i,k)-r(j,k)),2);
            }
            e_potential += 1/sqrt(r_12);
        }
    }
    return e_potential;
}


// LennardJones definitions
LennardJones::LennardJones(GeneralParams & gP) :
    Potential ( gP.number_particles, gP.dimension, "Lennard-Jones potential" ) {
        this->L = gP.L;
};

double LennardJones::get_pot_E(const mat& r) const {
    int i = 0, j = 0, k = 0;
    double e_potential = 0, r_12 = 0;
    for (i = 0; i < n_p-1; i++) {
        for (j = i+1; j < n_p; j++) {
            r_12 = 0;
            for (k = 0; k < dim; k++) {
                r_12 += pow( (r(i,k)-r(j,k)-round((r(i,k)-r(j,k))/L)*L), 2 );
            }
            e_potential += ( pow(r_12,-6) - pow(r_12,-3) );
        }
    }
    return 4*e_potential;
}
