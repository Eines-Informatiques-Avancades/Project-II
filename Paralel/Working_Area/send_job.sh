#!/bin/bash
# Run the program.
#$ -S /bin/bash
#$ -N procs55_parts13824
##$ -cwd
#$ -pe mpi 55
#$ -e name.err
time mpirun -np 55 program.exe

