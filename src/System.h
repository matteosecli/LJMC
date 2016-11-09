#include <vector>
#include <string>
#include <armadillo>

using namespace arma;

namespace QMC2
{

    class Potential;

    class System {
        protected:
            int n_p;
            int dim;
            std::vector<Potential*> potentials;

        public:
            System( int number_particles, int dimension );
            System();

            //! Method for adding a potential to the system.
            void add_potential(Potential* pot);

            //! Method for calculating the total potential energy.
            double get_potential_energy(const mat& r);
    };

    /* TODO: generalize with differnt types of systems, i.e. Fermions and Bosons.
     * One has to add
     *      class Fermions : public System { ... }
     */

}
