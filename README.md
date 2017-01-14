# LJMC


## Description

A Lennard-Jones Monte Carlo simulator written for the course "Understanding electron correlation with computer simulation" held at SISSA in the A.A. 2016-2017 by Sandro Sorella and Federico Becca.

Partially based on a VQMC simulator which I wrote for my BSc thesis ([https://github.com/matteosecli/QMC](https://github.com/matteosecli/QMC)).


## Installation

### Prerequisites
- `g++` (or `icpc`, not tested with `clang`)
- Armadillo ([http://arma.sourceforge.net/](http://arma.sourceforge.net/))
- LAPACK ([http://www.netlib.org/lapack/](http://www.netlib.org/lapack/))
- OpenBLAS or Intel MKL ([http://www.openblas.net/](http://www.openblas.net/) or [https://software.intel.com/en-us/intel-mkl](https://software.intel.com/en-us/intel-mkl))
- Libconfig ([http://www.hyperrealm.com/libconfig/](http://www.hyperrealm.com/libconfig/))
- `qmake` (just if you want to compile using the Makefile). Note that an installation of QT is not needed, you just need `qmake` (which doesn't depend on QT).

### Installation (with qmake)
Grab the latest version of **LJMC** from the GitHub repository:

	git clone https://github.com/matteosecli/LJMC && cd LJMC/src
	
Then, type

	qmake && make
	
If everything goes as it should, you now have an executable `LJMC` in your folder. You can run it via

	./LJMC outputfile
	
Refer to the [Usage](#usage) section for usage instructions.

### Installation (without qmake)
Grab the latest version of **LJMC** from the GitHub repository:

	git clone https://github.com/matteosecli/LJMC && cd LJMC/src
	
Then, type

	g++ -fPIC -c lib.cpp && g++ -o LJMC -std=c++11 -O2 lib.o main.cpp Potential.cpp System.cpp LennardJones.cpp pBar.cpp -lconfig++ -larmadillo -llapack -lblas

If everything goes as it should, you now have an executable `LJMC` in your folder. You can run it via

	./LJMC outputfile
	
Refer to the [Usage](#usage) section for usage instructions. Note that as soon as the program is modified, the `g++` compile command could change and this README could not be atply updated; so, I suggest to go the qmake/make way to be sure.


## Usage

### Passing the parameters
The program takes as the only input on command line the output file. You can run it via

	./LJMC outputfile
	
Then, a series of questions follow asking you the parameters of the simulation.

If you need to run a simulation without an interactive shell, you can pass the answers via e.g. a Bash script. Example:

	#!/bin/bash
	## WARNING: 'sh' does not work
	# Format:
	# -------
	# Number of particles = 100
	# Dimensionality = 3
	# Temperature = 2.0
	# Density = 0.5
	# Monte Carlo steps = 1000000
	# Step length = 0.32
	# Thermalization steps = 1000
	# Do you want to calculate the pressure? (1-y / 0-n) = 1
	# Do you want to calculate the g(r)? (1-y / 0-n) = 1
	
	./LJMC out_file  << PARAMETERS
	100
	3
	2.0
	0.5
	1000000
	0.32
	1000
	1
	1
	PARAMETERS
	
Saving the above code as `run.sh` and launching it via `./run.sh` will run the simulation without interactivly asking for user inputs.

In the future I'll try to improve the usability by giving the possibility of choosing between using an interactive shell and passing the required parameters via switches.


### Output files
Each run produces three different files: `outputfile`, `outputfile.cfg` and `outputfile.mat`.

- `outputfile` is the most important one. The first column contains the samples of the potential energy, while the second one contains the samples of the pressure (if the calculation of the pressure is required at input time).
- `outputfile.cfg` contains the informations needed to perfectly reproduce the system studied (e.g. number of paerticles, density, temperature, etc.) when continuing a calculation.
- `outputfile.mat` contains the positions of the particles at the end of the simulation. The file is in the `arma_binary` format (see [http://arma.sourceforge.net/docs.html#save_load_mat](http://arma.sourceforge.net/docs.html#save_load_mat)).

If you also asked for the calculation of the g(r), the program will produce an additional `g(r)_-_outputfile` file. The first column contains the variable `r`, while the second one contains the values of `g(r)`.


###Restart a calculation
This is actually more of a "*continue the calculation*" feature. It's useful if you realize that you've not produced enough Monte Carlo steps and you want to produce some more, without thermalizing again. The syntax in this case is

	./LJMC outputfile restart

where `outputfile` is the output file *produced by the calculation you want to continue*. The new samples will be appended at the end of the file. The program will look for a `outputfile.cfg` file and a `outputfile.mat` in order to load the previous calculation results, so please **don't change the names of the output files!**. If you still want to do so, change all the names accordingly.

Please note that the program *does not load* the pre-existing samples of potential energy and pressure; therefore, the results printed on screen at the end of the simulation only refers to the new samples. In order to analyze old and new samples together, just use the provided MATLAB scripts in the usual way.

Note also that the g(r), in this case, cannot be calculated. This is because I've still not instructed the program to save and reload the necessary data for the calculation of the g(r) in the restart mode; I hope I will find some time to implement that, eventually.