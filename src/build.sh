#g++ -fPIC -c lib.cpp
g++ -std=c++11 -fopenmp -O2 lib.o main.cpp Potential.cpp System.cpp LennardJones.cpp pBar.cpp -larmadillo -llapack -lblas

# For static:
# -static -static-libgcc -static-libstdc++

# For the cluster:
# I've already compiled armadillo. Unload gnu and load gcc 4.8
#g++ -O2 -fPIC -c lib.cpp 
#g++ -I ~/armadillo/usr/local/include/ -L ~/armadillo/usr/local/lib64/ -std=c++11 -fopenmp -O2 lib.o main.cpp Potential.cpp System.cpp LennardJones.cpp pBar.cpp -larmadillo
