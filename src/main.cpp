/*    CompPhys@SISSA - Lennard-Jones Project
 *          Matteo Seclì, Fall 2016
 *          <secli.matteo@gmail.com>
 *
 *
 * This is the first project for the course "Understanding
 * electron correlation with computer simulation" held at
 * SISSA in the AA 2016-2017 by Sandro Sorella and Federico
 * Becca.
 *
 * It's a MC simulator for Lennard-Jones fluids, which by
 * design choice works with any number of particles and in
 * any dimension (from 1D to 3D, even if I've not tested it
 * in 1D). The case of interest in this project is a fluid
 * of 100 particles in a cubic box with PBC's, such that
 * the density is 0.5. That's why the density is fixed in
 * the 'initialize()' function; you can easily ask for a
 * user input if you need that.
 *
 * The code recycles parts of a previous project for my BSc
 * thesis, developed at the University of Oslo under the
 * supervision of Morten Hjorth-Jensen. It consisted of a
 * VMC simulator, that's why there are bits af variational
 * traces in the code. I've tried to keep the original
 * object structure; in 'Potential.cpp' there are now some
 * additional definitions for a LJ potential that are not
 * acually used in the final code, since it's specifically
 * optimized for the LJ with OPM's. Anyway, they are there
 * for you if you want e.g. to move all the particles
 * together or if you want to use a different potential.
 * The description of the previous project is left below
 * for reference.
 *
 * TODO: (in descending order of importance)
 *  - Optimize memory usage by using a sparse matrix for
 *    the potential.
 *  - Make the specific MC functions more general and
 *    transfer them in a different file.
 *  - Give the ability to use different boundary
 *    conditions, e.g. OBC. This seems a job tailored
 *    to the object structure nature of the project.
 *
 *
 * The code is hosted at:
 *  <https://github.com/matteosecli/LJMC>
 *
 * Refer to the LICENSE file for license information.
 *
 */


/* -------------------- OLD HEADER -------------------- */
/*          FYS3150@UiO - QMC Project
 *          Matteo Seclì, Autumn 2014
 *          <secli.matteo@gmail.com>
 *
 * This program performs a Variational Monte Carlo
 * calculation for Quantum Dots. It works both in 2
 * and 3 dimensions, but you have to use a number
 * of electrons equal to the degeneracy of the energy
 * level you want to consider.
 * The program is structured in an object-oriented way,
 * in order to perform further generalization after
 * this project. As an example, you can easily remove
 * electron-electron interaction (just comment one line)
 * or change the program to perform calculations over
 * an atom instead of an harmonic oscillator (again,
 * just one line to change).
 *
 * The main object structure is thanks to the work of
 * Jørgen Høgberget, from whom I've borrowed many lines
 * of code. He's done a great job with his project.
 * You can check his work here:
 *  <https://github.com/jorgehog/QMC2>
 * and here:
 *  <http://urn.nb.no/URN:NBN:no-38645>
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include "LennardJones.h"
#include <armadillo>
#include "Potential.h"
#include "System.h"
#include "Structs.h"
#include <time.h>
#include <vector>
#include "pBar.h"

using namespace  std;
using namespace arma;
using namespace QMC2;
using namespace pBarNamespace;

/* Output file as global variable */
ofstream ofile;

/* The step length and its squared inverse for the second derivative */
#define h 0.001
#define h2 1000000



/*--- Declaraton of functions ---*/

/* Function to read in data from screen, note call by reference */
void initialize(GeneralParams&, VMCparams&) ;

/* The MC sampling for the variational Monte Carlo */
void mc_sampling(GeneralParams&, VMCparams&, double&, double&, System&);

/* The Metropolis algo performed in a brute force way */
void metropolis_bf(GeneralParams&, VMCparams&, System&, mat&, mat&,
                   double&, double&, int&);



/*--- Begin of main program ---*/
int main(int argc, char* argv[])
{
    /* Start timing */
    //boost::timer t;

    /* Initialize variables */
    char *outfilename;
    double cumulative_e = 0, cumulative_e2 = 0;

    /* Read in output file, abort if there are too few
     * command-line arguments.
     */
    if( argc <= 1 ){
        cout << "Bad Usage: " << argv[0] << " read also output file on same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }
    ofile.open(outfilename);

    /* Initialize structs */
    struct GeneralParams gP;
    struct VariationalParams vP;
    struct VMCparams vmcParams;

    /* Initialize system parameters */
    initialize(gP, vmcParams);

    /* Build the required system */
    System MySystem;
    MySystem = System(gP.number_particles, gP.dimension);
    MySystem.add_potential(new LennardJones(gP));
    /* To be generalized into (e.g.) MySystem = new Fermions(gP).
     * After generalization, substitute with
     *      MySystem->add_potential(new Harmonic_osc(gP))
     */

    /* Do the MC sampling */
    mc_sampling(gP, vmcParams, cumulative_e, cumulative_e2, MySystem);

    /* Close the output file */
    ofile.close();

    /* Print out elapsed time */
    //double time = t.elapsed();
    //cout <<"Elapsed time:\t\t" << time <<"s." <<endl;

    return 0;
}


/* Monte Carlo sampling with the Metropolis algorithm */
void mc_sampling(GeneralParams& gP, VMCparams& vmcParams,
                 double& cumulative_e, double& cumulative_e2, System& InputSystem)
{
    int accept;
    double energy, energy2;
    int npartperL;
    double Lspacing;
    mat r_old, r_new;

    /* Allocate matrices which contain the position of the particles */
    r_old = zeros<mat>(gP.number_particles, gP.dimension);
    r_new = zeros<mat>(gP.number_particles, gP.dimension);

    /* Initialize the matrices such that we have no divergencies */
    npartperL = (int) ceil(pow(gP.number_particles, 1.0/((double) gP.dimension)));
    Lspacing = gP.L/((double) npartperL);
    /* Leave the first particle at r=0 */
    for ( int i = 1; i < gP.number_particles; i++ )
    {
        /* x_i */
        r_old(i,0) = (i%npartperL)*Lspacing;
        /* y_i and z_i */
        for ( int k = 1; k < gP.dimension; k++ )
        {
            r_old(i,k) = ((int) i / ((int) pow(npartperL,k)))*Lspacing;
        }
    }

    r_new = r_old;

    /* Initialize the V_old matrix */
    mat V_old = zeros<mat>(gP.number_particles,gP.number_particles);
    LJPotMatrix(gP, V_old, r_old);
    cout << "Initial energy = " << LJPotValue(gP, V_old)+LJPotTail(gP) << endl;

    /* Initialization of the energies */
    energy = energy2 = 0; accept =0;

    /* Perform Metropolis sampling */
    metropolis_bf(gP, vmcParams, InputSystem, r_old, r_new,
                  energy, energy2, accept);

    /* Update the energy average and its square */
    cumulative_e = energy/vmcParams.number_cycles;
    cumulative_e2 = energy2/vmcParams.number_cycles;

    /* Calculate some optimistic error and the Cv */
    double variance, error, AcceptPercent, Cv;
    variance = cumulative_e2-pow(cumulative_e,2);
    error=sqrt(variance/vmcParams.number_cycles);
    Cv = variance/pow(gP.temp,2);
    AcceptPercent = ((double) accept)/((double) vmcParams.number_cycles )*100.0;
    cout << "V = " << setprecision(8) << cumulative_e << " "
         << "+/- " << setprecision(8) << error << "\t"
         << "Cv = " << setprecision(8) << Cv << "\t"
         << "Accepted steps = " << setprecision(2) << AcceptPercent << "%" << endl;

    /* Free memory */
    r_old.reset(); // free memory
    r_new.reset(); // free memory
}   /* END of mc_sampling function */


/* The Metropolis brute-force algo */
void metropolis_bf(GeneralParams& gP, VMCparams& vmcParams,
                   System& InputSystem, mat& r_old, mat& r_new,
                   double& energy, double& energy2, int& accept)
{
    /* Initialize some needed variables */
    int moving_particle;
    long idum = -1;
    double delta_e = 0;
    double mratio = 0;
    double BoltzmannBeta = 1.0/gP.temp;
    mat V_old, V_new;
    double PotDiff = 0;
    double PotTailCorrection = 0;
    int PrintRate = 1000;
    double ProgressIncrement = 100.0*((double) PrintRate)/((double) vmcParams.number_cycles+vmcParams.thermalization);
    pBar progressbar;

    /* Matrices for the potential. Maybe it's a waste of space, but we have
     * plenty of it and I'm lazy. I could optimize for memory by using
     * the 'sp_mat' typedef from Armadillo.
     */
    V_old = zeros<mat>(gP.number_particles,gP.number_particles);
    V_new = zeros<mat>(gP.number_particles,gP.number_particles);
    LJPotMatrix(gP, V_old, r_old);
    V_new = V_old;

    /* Calculate the tail correction to the potential,
     * which is always the same. */
    PotTailCorrection = LJPotTail(gP);

    /* Loop over the Monte Carlo cycles */
    for (int cycles = 1; cycles <= vmcParams.number_cycles+vmcParams.thermalization; cycles++){

        /* Propose a new position, move only
         * one particle at a time!
         */
        moving_particle = round(ran1(&idum)*(gP.number_particles-1));
        for ( int j=0; j < gP.dimension; j++ ) {
            r_new(moving_particle,j) = r_old(moving_particle,j)+vmcParams.step_length*(ran1(&idum)-0.5);
        }

        /* Calculate the acceptance ratio by the OPM optimization.
         * Fistly, one puts V_new equal to V_old, because then V_new
         * is gonna be modified. Then, delta_e is assigned the value
         * V_old has. Finally, LJDiffValue updates V_new to its actual
         * value and returns the difference (V_new-V_old). The Metropolis
         * acceptance ratio is then calculated.
         */
        V_new = V_old;
        delta_e = LJPotValue(gP, V_old); // So far V_new is V_old
        PotDiff = LJDiffValue(gP, V_new, r_new, moving_particle); // This function modifies V_new into the actual V_new
        mratio = exp(-PotDiff*BoltzmannBeta);

        /* Metropolis test */ //< or <=?
        if ( (mratio < 1) ? (ran1(&idum) < mratio) : 1 ) {
            for (int  j=0; j < gP.dimension; j++) {
                r_old(moving_particle,j)=r_new(moving_particle,j);
            }
            delta_e += PotDiff;
            V_old = V_new;
            if ( cycles > vmcParams.thermalization ) { accept = accept+1; }
        }

        /* Compute the energy */
        if ( cycles > vmcParams.thermalization ) {
            /* Add the tail correction to the potential */
            delta_e += PotTailCorrection;
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << setprecision(8) << delta_e << endl;
            energy += delta_e;
            energy2 += pow(delta_e,2);
        }

        /* Print some nice visual feedback every PrintRate points
         */
        if ( cycles%PrintRate == 0 ) {
            /* Update pBar by giving to it the current percentage */
            progressbar.update(ProgressIncrement);
            /* Print pBar */
            progressbar.print();
        }

    }   /* END of loop over MC trials */

    /* Print final percentage of the progress bar
     */
    progressbar.update(ProgressIncrement);
    progressbar.print();
    cout << endl;

}   /* END of metropolis_bf function */


/* Initialize the parameters needed for the simulation */
void initialize(GeneralParams & gP, VMCparams & vmcParams)
{
    cout << "Number of particles = ";
    cin >> gP.number_particles;

    cout << "Dimensionality = ";
    cin >> gP.dimension;

    cout << "Temperature = ";
    cin >> gP.temp;

    /* Set the density */
    cout << "Density = ";
    cin >> gP.density;

    /* Be sure to deactivate the parallel capabilities */
    gP.num_threads = 1;

    cout << "# Monte Carlo steps = ";
    cin >> vmcParams.number_cycles;

    cout << "# Step length = ";
    cin >> vmcParams.step_length;

    cout << "# Thermalization steps = ";
    cin >> vmcParams.thermalization;

    /* Turn off the variational capabilities */
    vmcParams.max_variations = 0;

    /* Set the box length */
    gP.L = pow(((double) gP.number_particles)/gP.density, 1.0/((double) gP.dimension));

    /* Set the square of the cutoff */
    gP.potcutoff2 = gP.L*gP.L/4.0;

    // DEBUG
    //cout << gP.potcutoff2 << endl;

}  // END of function initialize
