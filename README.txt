This Skyrme Hartree-Fock code is modularized into 5 files:
pnInfinite.py (Main file)
parameters.py
npSubroutines.py
states.py
fcoul2py.f90

pnInfinite.py is the main file and contains the basic algorithm for solving HF equation
parameters.py sets constants and other parameters used in the rest of the code; it also defines command line arguements
states.py defines the State object used for holding eigenfunctions of the matrix
npSubroutines.py contains all helper functions EXCEPT those relating to the Coulomb potential

fcoul2py.f90 is a fortran subroutine designed to combile into a single python function for solving the Coulomb potential.
This was done to optimize speed and memory and allow parallelization through OpenMP.

To run the algorithm, there are 2 steps

First: compile fortran subroutine file into python function, with openmp, using the command
f2py -c -m fcoul2py fcoul2py.f90 --f90flags="-fopenmp" -lgomp
Note: you will need module scipy-stack loaded on SHARCNET

Second: run the python interpreter on the main file with the following command line arguements
python pnInfinite.py nMult cubeMult YP NXNYMAX NQ rootpath

where
nMult * (cubeMult^3) = total number of particles
yp is fraction of particles that are protons, given as integer i.e. 50 is 50%
NXNYMAX is max nx ny space to explore
NQ is number of periods in external potential (vq is set in parameters.py)
rootpath directory where you set MOST of the output to go; data written to files, but not stout

The output of this program is several files relating to density written to rootpath, 
as well as line printed to stout that displays the final E/A value.
This can be seen towards the end of the main file.
