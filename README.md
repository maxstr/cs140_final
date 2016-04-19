CS140 Final Project
====
We're approximating a heat equation for a given heat source term. Here it's sin(x) + sin(y).


Usage
====
To use:

set p, rootp, m, n in twoDOutput.cpp. M must equal N and be multiples of rootP

compile as mpicc -o binaryName twoDOutput.cpp

mpirun -machinefile $PBS_NODEFILE -np p binaryName
