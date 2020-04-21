Build instruction
===============================
In serial
module load gcc 

gcc K_means_iter.c -lm
===============================
In parallel the files for
source load--------------------loads the necesaary modules.
source compile-----------------compiles the code and creates executable.

Assignment 4 example
=================================
module load gcc/7.4.0/1 spectrum-mpi
mpic++ <source file>
mpirun -np 10 <executable>
