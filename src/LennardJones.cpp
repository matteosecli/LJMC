#include "LennardJones.h"

using namespace arma;
using namespace QMC2;

/* Returns the potential matrix V calculated at r */
void LJPotMatrix(GeneralParams & gP, mat& V, const mat& r) {
    V.zeros();
    int i = 0, j = 0, k = 0;
    double r_12 = 0;
    double temp = 0;
    for (i = 0; i < gP.number_particles -1; i++) {
        for (j = i+1; j < gP.number_particles; j++) {
            r_12 = 0;
            for (k = 0; k < gP.dimension; k++) {
                r_12 += pow( (r(i,k)-r(j,k)-round((r(i,k)-r(j,k))/gP.L)*gP.L), 2 );
            }
            /* Do a cutoff at L/2 */
            if ( r_12 < gP.potcutoff2 ) {
                temp = pow(r_12,-3);
                V(i,j) = 4*( pow(temp,2) - temp );
            } else {
                V(i,j) = 0;
            }
        }
    }
}


/* Returns the value of the potential matrix V
 * by exploiting the fact that it's upper triangular
 */
double LJPotValue(GeneralParams & gP, const mat& V) {
    int i = 0, j = 0;
    double V_value = 0;
    for (i = 0; i < gP.number_particles -1; i++) {
        for (j = i+1; j < gP.number_particles; j++) {
            V_value += V(i,j);
        }
    }
    return V_value;
}


/* Returns the difference of the potential matrix V after updating
 * the elements relative to the moving particle. WARNING! mov_part
 * is already the INDEX of the moving particle (i.e. from 0 to n_p-1),
 * so no need to subtract 1.
 */
double LJDiffValue(GeneralParams & gP, mat& V, const mat& r, const int& mov_part) {
    int i = 0, j = 0, k = 0;
    double V_buffer = 0, DeltaV = 0, r_12 = 0;
    double temp = 0;
    for (i = 0; i < mov_part; i++) {
        r_12 = 0;
        for (k = 0; k < gP.dimension; k++) {
            r_12 += pow( (r(i,k)-r(mov_part,k)-round((r(i,k)-r(mov_part,k))/gP.L)*gP.L), 2 );
        }
        /* Do a cutoff at L/2 */
        if ( r_12 < gP.potcutoff2 ) {
            temp = pow(r_12,-3);
            V_buffer = 4*( pow(temp,2) - temp );
        } else {
            V_buffer = 0;
        }
        DeltaV += V_buffer - V(i,mov_part);
        V(i,mov_part) = V_buffer;
    }
    for (j = mov_part+1; j < gP.number_particles; j++) {
        r_12 = 0;
        for (k = 0; k < gP.dimension; k++) {
            r_12 += pow( (r(j,k)-r(mov_part,k)-round((r(j,k)-r(mov_part,k))/gP.L)*gP.L), 2 );
        }
        /* Do a cutoff at L/2 */
        if ( r_12 < gP.potcutoff2 ) {
            temp = pow(r_12,-3);
            V_buffer = 4*( pow(temp,2) - temp );
        } else {
            V_buffer = 0;
        }
        DeltaV += V_buffer - V(mov_part,j);
        V(mov_part,j) = V_buffer;
    }
    return DeltaV;
}


/* Calculates the tail correction to the potential */
double LJPotTail(GeneralParams & gP) {
    return 8.0/3.0*datum::pi*gP.number_particles*gP.density*(1/3*pow(gP.potcutoff2,-4.5)-pow(gP.potcutoff2,-1.5));
}
