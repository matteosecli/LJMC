#ifndef POTENTIAL_H
#define POTENTIAL_H

using namespace arma;

namespace QMC2
{

    struct GeneralParams;

    class Potential {
        protected:
            int n_p;
            int dim;
            std::string name;
        public:
            Potential ( int number_particles , int dimension, std::string name = "Potential" );
            Potential();
            // Pure virtual function
            virtual double get_pot_E ( const mat& r ) const = 0;
            std::string get_name(){return name;}
    };

    // Harmonic oscillator potential
    class Harmonic_osc : public Potential {
        protected:
           double w;
        public:
            Harmonic_osc (GeneralParams &) ;    // Here, 'frequency' is the 'omega'
            double get_pot_E (const mat &r ) const ;
    };

    // Coulomb potential in an atom
    class Coulomb : public Potential {
        protected :
            int Z;
        public :
            Coulomb (GeneralParams &) ;   // Remember that charge is in natural units: e=1
            double get_pot_E ( const mat& r ) const ;
    };

    // Repulsion potential between electrons
    class eRepulsion : public Potential {
        public :
            eRepulsion (GeneralParams &) ;
            double get_pot_E ( const mat& r ) const ;
    };

    // Lennard-Jones potential
    class LennardJones : public Potential {
        protected :
            double L;
        public :
            LennardJones (GeneralParams &) ;
            double get_pot_E (const mat& r) const ;
    };
}

#endif //POTENTIAL_H
