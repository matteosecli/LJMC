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

	g++ -fPIC -c lib.cpp && g++ -o LJMC -std=c++11 -O2 lib.o main.cpp Potential.cpp System.cpp LennardJones.cpp pBar.cpp -larmadillo -llapack -lblas

If everything goes as it should, you now have an executable `LJMC` in your folder. You can run it via

	./LJMC outputfile
	
Refer to the [Usage](#usage) section for usage instructions. Note that as soon as the program is modified, the `g++` compile command could change and this README could not be atply updated; so, I suggest to go the qmake/make way to be sure.


## Usage

The program takes as the only input on command line the output file. You can run it via

	./LJMC outputfile
	
Then, a series of questions follow asking you the parameters of the simulation. The density is fixed at 0.5, therefore it is not asked.

If you need to run a simulation without an interactive shell, you can pass the answers via e.g. a Bash script. Example:

	# Format:
	# -------
	# Number of particles = 100
	# Dimensionality = 3
	# Temperature = 2.0
	# Monte Carlo steps = 1000000
	# Step length = 0.32
	# Thermalization steps = 1000
	
	./LJMC out_file  << PARAMETERS
	100
	3
	2.0
	1000000
	0.32
	1000
	PARAMETERS
	
Saving the above code as `run.sh` and launching it via `./run.sh` will run the simulation without interactivly asking for user inputs.

In the future I'll try to improve the usability by giving the possibility of choosing between using an interactive shell and passing the required parameters via switches.