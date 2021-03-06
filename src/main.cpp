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
#include <string>
#include "lib.h"
#include "LennardJones.h"
#include <armadillo>
#include "Potential.h"
#include "System.h"
#include "Structs.h"
#include <time.h>
#include <vector>
#include "pBar.h"
#include <libconfig.h++>

using namespace  std;
using namespace arma;
using namespace libconfig;
using namespace QMC2;
using namespace pBarNamespace;

/* Output file as global variable */
ofstream ofile;
char *outfilename;

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
                   double&, double&, double&, double&, int&);



/*--- Begin of main program ---*/
int main(int argc, char* argv[])
{
    /* Start timing */
    //boost::timer t;

    /* Initialize variables */
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

    /* Initialize structs */
    struct GeneralParams gP;
    struct VariationalParams vP;
    struct VMCparams vmcParams;

    /* Determine whether the calculations is a fresh one
     * or a restarted one, and open the output file accordingly
     */
    if ( argc == 2 ) {
        ofile.open(outfilename);
    } else if ( argc == 3 && string(argv[2]) == "restart" ) {
        gP.restart = true;
        ofile.open(outfilename, std::ofstream::in | std::ofstream::app);
        cout << "Calculation is restarting from where it left..." << endl;
    } else {
        cout << "Bad Usage: unknown option '" << argv[2] << "'. ";
        cout << "Use 'restart' (without quotes) if you want to restart the calculation" << endl;
    }

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
    double pressure, pressure2;
    int npartperL;
    double Lspacing;
    mat r_old, r_new;

    /* Allocate matrices which contain the position of the particles */
    r_old = zeros<mat>(gP.number_particles, gP.dimension);
    r_new = zeros<mat>(gP.number_particles, gP.dimension);

    /* Initialize the matrices such that we have no divergencies */
    if ( gP.restart ) {
        r_old.load(string(outfilename)+".mat");
    } else {
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
    }

    r_new = r_old;

    /* Initialize the V_old matrix */
    mat V_old = zeros<mat>(gP.number_particles,gP.number_particles);
    LJPotMatrix(gP, V_old, r_old);

    /* Initialization of the energies */
    energy = energy2 = 0; accept =0;

    /* Initialization of pressures */
    pressure = pressure2 = 0;

    /* Perform Metropolis sampling */
    metropolis_bf(gP, vmcParams, InputSystem, r_old, r_new,
                  energy, energy2, pressure, pressure2, accept);

    /* Update the energy average and its square */
    cumulative_e = energy/vmcParams.number_cycles;
    cumulative_e2 = energy2/vmcParams.number_cycles;

    /* Normalize the pressure */
    pressure = pressure/vmcParams.number_cycles;
    pressure2 = pressure2/vmcParams.number_cycles;

    /* Calculate some optimistic error and the Cv */
    double variance, error, AcceptPercent, Cv;
    variance = cumulative_e2-pow(cumulative_e,2);
    error=sqrt(variance/vmcParams.number_cycles);
    Cv = variance/pow(gP.temp,2);
    AcceptPercent = ((double) accept)/((double) vmcParams.number_cycles )*100.0;
    cout << "V/N = " << setprecision(8) << cumulative_e/gP.number_particles << " "
         << "+/- " << setprecision(8) << error/gP.number_particles << "\t"
         << "Cv/N = " << setprecision(8) << Cv/gP.number_particles << "\t";
         if (vmcParams.DoPressure) {
            cout << "P = " << setprecision(8) << pressure << " "
            << "+/- " << setprecision(8) << sqrt((pressure2-pow(pressure,2))/vmcParams.number_cycles) << endl;
         }
         cout << "Accepted steps = " << setprecision(2) << AcceptPercent << "%" << endl;

    /* Save the positions of the particles */
    r_new.save(string(outfilename)+".mat");

    /* Free memory */
    r_old.reset(); // free memory
    r_new.reset(); // free memory
}   /* END of mc_sampling function */


/* The Metropolis brute-force algo */
void metropolis_bf(GeneralParams& gP, VMCparams& vmcParams,
                   System& InputSystem, mat& r_old, mat& r_new,
                   double& energy, double& energy2,
                   double& pressure, double& pressure2,
                   int& accept)
{
    /* Initialize some needed variables */
    int moving_particle;
    long idum = -1;
    double delta_e = 0;
    double delta_p = 0;
    double mratio = 0;
    double BoltzmannBeta = 1.0/gP.temp;
    mat V_old, V_new;
    mat W_old;
    double PotDiff = 0;
    double PotTailCorrection = 0;
    double VirTailCorrection = 0;
    double PressExtValue = 0;
    int PrintRate = 1000;
    double ProgressIncrement = 100.0*((double) PrintRate)/((double) vmcParams.number_cycles+vmcParams.thermalization);
    pBar progressbar;

    /* Variables for g(r) */
    int gRate = 10000;
    int NgTrials = 0;
    int NgBins = 200;
    int iBin = 0;
    double gBinWidth = gP.L/(2.0*((double)NgBins));
    double gDistance = 0;
    vec gr = zeros<vec>(NgBins);

    /* Matrices for the potential. Maybe it's a waste of space, but we have
     * plenty of it and I'm lazy. I could optimize for memory by using
     * the 'sp_mat' typedef from Armadillo.
     */
    V_old = zeros<mat>(gP.number_particles,gP.number_particles);
    V_new = zeros<mat>(gP.number_particles,gP.number_particles);
    LJPotMatrix(gP, V_old, r_old);
    V_new = V_old;
    /* Do the same for the virial W */
    W_old = zeros<mat>(gP.number_particles,gP.number_particles);
    LJVirMatrix(gP, W_old, r_old);

    /* Calculate the tail correction to the potential,
     * which is always the same. */
    PotTailCorrection = LJPotTail(gP);
    /* Do the same for the virial W */
    VirTailCorrection = LJVirTail(gP);
    PressExtValue = gP.density*gP.temp;

    if ( vmcParams.tailcorr_active ) {
        cout << "Using tail corrections: \t"
        << "ΔV/N = " << PotTailCorrection/gP.number_particles << "\t"
        << "and" << "\t"
        << "ΔP = " << VirTailCorrection << endl;
    } else {
        cout << "Not using tail corrections." << endl;
    }

    /* Calculate the initial energy */
    delta_e = LJPotValue(gP, V_old) + ( vmcParams.tailcorr_active ? PotTailCorrection : 0.0 );

    /* Calculate the initial pressure */
    delta_p = PressExtValue + LJVirValueOverV(gP, W_old) + ( vmcParams.tailcorr_active ? VirTailCorrection : 0.0 );

    /* Print some initial information */
    cout << "Initial energy per particle = "
         << delta_e/gP.number_particles << endl;
    if ( vmcParams.tailcorr_active ) {
        cout << "Initial pressure = "
        << delta_p << endl;
    }

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
        /* SF -- Slightly safer (less error propagation) but slower method. */
        //delta_e = LJPotValue(gP, V_old); // So far V_new is V_old
        PotDiff = LJDiffValue(gP, V_new, r_new, moving_particle); // This function modifies V_new into the actual V_new
        mratio = exp(-PotDiff*BoltzmannBeta);

        /* Metropolis test */
        if ( (mratio < 1) ? (ran1(&idum) < mratio) : 1 ) {
            /* Update the virial if asked */
            if (vmcParams.DoPressure) {
                delta_p = PressExtValue + LJUpdateVirValueOverV(gP, W_old, r_new, moving_particle) + ( vmcParams.tailcorr_active ? VirTailCorrection : 0.0 );
            }
            /* Accept the move */
            for (int  j=0; j < gP.dimension; j++) {
                r_old(moving_particle,j)=r_new(moving_particle,j);
            }
            /* Update the energy */
            delta_e += PotDiff;
            V_old = V_new;
            /* Update the acceptance countings */
            if ( cycles > vmcParams.thermalization ) { accept = accept+1; }
        } else {
            /* Reject the move */
            for (int  j=0; j < gP.dimension; j++) {
                r_new(moving_particle,j)=r_old(moving_particle,j);
            }
        }

        /* Compute the energy */
        if ( cycles > vmcParams.thermalization ) {
            /* Add the tail correction to the potential. Only needed
             * when using SF (safer but slower method). */
            //delta_e += PotTailCorrection;

            /* Fill the output file and the cumulants */
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << setprecision(8) << delta_e;
            if (vmcParams.DoPressure) {
                ofile << "\t" << setprecision(8) << delta_p;
                pressure += delta_p;
                pressure2 += pow(delta_p,2);
            }
            ofile << endl;
            energy += delta_e;
            energy2 += pow(delta_e,2);
        }

        /* Calculate g(r) every gRate points */
        if ( cycles > vmcParams.thermalization && cycles%gRate == 0 && vmcParams.DoGr ) {
            for ( int i = 0; i < gP.number_particles - 1; i++ ) {
                for ( int j = i+1; j < gP.number_particles; j++ ) {
                    gDistance = 0;
                    for ( int k = 0; k < gP.dimension; k++ ) {
                        gDistance += pow( (r_new(i,k)-r_new(j,k)-round((r_new(i,k)-r_new(j,k))/gP.L)*gP.L), 2 );
                    }
                    if ( gDistance < gP.potcutoff2 ) {
                        gDistance = sqrt(gDistance);
                        iBin = floor(gDistance/gBinWidth);
                        if ( iBin == NgBins ) iBin -= 1;
                        gr(iBin) += 2;
                    }
                }
            }
            NgTrials += 1;
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

    /* Renormalize g(r) and print it */
    if ( vmcParams.DoGr ) {
        ofstream grofile;
        grofile.open("g(r)_-_"+string(outfilename));
        for ( int i = 0; i < NgBins; i++ ) {
            gr(i) = 3.0*gr(i)/(4.0*datum::pi*((double)NgTrials)*gP.density*((double)gP.number_particles)*pow(gBinWidth,3.0)*((double)(pow(i+2,3.0)-pow(i+1,3))));
            grofile << setprecision(8) << (i+1)*gBinWidth << "\t" << gr(i) << endl;
        }
        grofile.close();
    }

    /* Print final percentage of the progress bar
     */
    progressbar.update(ProgressIncrement);
    progressbar.print();
    cout << endl;

}   /* END of metropolis_bf function */


/* Initialize the parameters needed for the simulation */
void initialize(GeneralParams & gP, VMCparams & vmcParams)
{
    int press_active_resp = 0;
    int gr_active_resp = 0;

    /* Configuration file filename */
    string cfgoutfilename = string(outfilename) + ".cfg";

    /* Configuration file object */
    Config conf;

    if ( !gP.restart ) {
        cout << "Number of particles = ";
        cin >> gP.number_particles;

        cout << "Dimensionality = ";
        cin >> gP.dimension;

        cout << "Temperature = ";
        cin >> gP.temp;

        /* Set the density */
        cout << "Density = ";
        cin >> gP.density;

        cout << "# Monte Carlo steps = ";
        cin >> vmcParams.number_cycles;

        cout << "# Step length = ";
        cin >> vmcParams.step_length;

        cout << "# Thermalization steps = ";
        cin >> vmcParams.thermalization;

        cout << "Do you want to calculate the pressure? (1-y / 0-n) = ";
        cin >> press_active_resp;
        if (press_active_resp != 0) vmcParams.DoPressure = true;

        cout << "Do you want to calculate the g(r)? (1-y / 0-n) = ";
        cin >> gr_active_resp;
        if (gr_active_resp != 0) vmcParams.DoGr = true;

        /* Build the configuration */
        Setting &root = conf.getRoot();
        root.add("number_particles", Setting::TypeInt) = gP.number_particles;
        root.add("dimension", Setting::TypeInt) = gP.dimension;
        root.add("temp", Setting::TypeFloat) = gP.temp;
        root.add("density", Setting::TypeFloat) = gP.density;
        root.add("number_cycles", Setting::TypeInt) = vmcParams.number_cycles;
        root.add("step_length", Setting::TypeFloat) = vmcParams.step_length;
        root.add("DoPressure",Setting::TypeBoolean) = vmcParams.DoPressure;

    } else {
        cout << "Loading configuration from config file..." << endl;

        /* Open the configuration file */
        conf.readFile(cfgoutfilename.c_str());

        /* Load previous configuration */
        int oldMCSteps;
        gP.number_particles = conf.lookup("number_particles");
        gP.dimension = conf.lookup("dimension");
        gP.temp = conf.lookup("temp");
        gP.density = conf.lookup("density");
        oldMCSteps = conf.lookup("number_cycles");
        vmcParams.step_length = conf.lookup("step_length");
        vmcParams.DoPressure = conf.lookup("DoPressure");

        /* Print info about the current configuration */
        cout << "You are using the following configuration:" << endl;
        cout << "Number of particles = " << gP.number_particles << endl;
        cout << "Dimensionality = " << gP.dimension << endl;
        cout << "Temperature = " << gP.temp << endl;
        cout << "Density = " << gP.density << endl;
        cout << "# Step length = " << vmcParams.step_length << endl;
        cout << "Calculation of the pressure = " << vmcParams.DoPressure << endl;

        /* Ask for additional MC steps */
        cout << "You've already done " << oldMCSteps << " Monte Carlo steps." << endl;
        cout << "How many additional Monte Carlo steps would you like to do?" << endl;
        cout << "# Additional Monte Carlo steps = ";
        cin >> vmcParams.number_cycles;

        /* Calculation of g(r) in restart mode is currently unavailable */
        vmcParams.DoGr = false;

        /* We don't need anymore to thermalize */
        vmcParams.thermalization = 0;

        /* Build the configuration. No need to rewrite old conf. */
        Setting &root = conf.getRoot();
        root["number_cycles"] = vmcParams.number_cycles + oldMCSteps;

    }

    /* Be sure to deactivate the parallel capabilities */
    gP.num_threads = 1;

    /* Turn off the variational capabilities */
    vmcParams.max_variations = 0;

    /* Set the box length */
    gP.L = pow(((double) gP.number_particles)/gP.density, 1.0/((double) gP.dimension));

    /* Set the square of the cutoff */
    gP.potcutoff2 = gP.L*gP.L/4.0;

    /* Save the configuration */
    conf.writeFile(cfgoutfilename.c_str());

}  // END of function initialize
