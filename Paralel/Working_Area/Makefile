SHELL := /bin/bash # necessary to use "time" shell command
# Parallel version of the Makefile for the Molecular Dynamics simulation

# MPI Compiler

mpi_comp = mpif90

# MPI run

mpi_exe = mpirun

# Some path and object variables are defined

OTHERDIR_MAIN=../Code

OTHERDIR_MOD=../Code/Modules

OTHERDIR_RESULTS=../Results/Data

OTHERDIR_FIG=../Results/Figures

MOD_ORDENADOS=../Code/Modules/Parallel_module.f90 ../Code/Modules/Constants_module.f90 ../Code/Modules/Pbc_module.f90 ../Code/Modules/Initial_state_module.f90 ../Code/Modules/Forces_module.f90 ../Code/Modules/Integrator_module.f90 ../Code/Modules/Gofr_module.f90 ../Code/Main.f90

#----------------------------------------------------------------------------------------------------------------------------------#

## Information about how to use this makefile.
##
## Option 1: After entering the parameters in parameters.txt, you can use make all to realize the whole process: Simulation + 
##           Statistical analysis + Figures.
##           If you use this method, the data will automatically be saved in the folder Results/Data and the figures in the folder 
##           Results/Figures. Additionally, the unnecessary files will be removed.
##
## Option 2: After entering the parameters in parameters.txt, you can use the commands that are explained below to realize the process 
##           step by step (less practical in my opinion).
##           If you want to do this, you must keep in mind that you only can analyse files that are in this folder, so do not use the 
##           move command until you have analysed the files.
##           Additionally, statistical analysis programs work with the files that are specified in parameters.txt, so make sure the 
##           file names are correct.
##
## Send Job: If you want to submit a job to a remote cluster, you can use this option (be sure the cluster has the MPI module). 
##           After entering the parameters in parameters.txt and properly modify the script send_job.sh (number of processors 
##           and job name basically), write "Make send_job". This command will send the job (only the simulation part).
##
#----------------------------------------------------------------------------------------------------------------------------------#

## Commands

## all: compile fortran code, run the code, realize the statistical analysis and the figures, move the results to their folder and remove the useless files.

all: compile exe statistics move clean

## run: compiles the fortran code and runs the program.

run: compile exe

## exe: runs the program (The number of processors to use is requested.)

exe:
	source execution.sh

## send_job: compiles the fortran code and send a job. (The number of processors to be used must be introduced in "send_job.sh")

send_job: compile exe_cluster

## exe_cluster: send a job.

exe_cluster:
	qsub send_job.sh

## compile: compiles the fortran code

compile: $(MOD_ORDENADOS)
	@$(mpi_comp) -g $(MOD_ORDENADOS) -o program.exe


## statistics: Statistical analysis is performed and the figures are generated. You will be asked some questions.
.PHONY: statistics
statistics:
	@python $(OTHERDIR_MAIN)/generate_figures.py parameters.txt
	@python $(OTHERDIR_MAIN)/Block_avarage.py parameters.txt

## move: Results files (extensions .dat and .png) are moved to Results folders (Data and Figures).
.PHONY: move	
move: 
	@mv *.dat $(OTHERDIR_RESULTS)
	@mv *.png $(OTHERDIR_FIG)
	@mv *.xyz $(OTHERDIR_RESULTS)

## clean: Objects, executables and modules are removed.
.PHONY: clean
clean:
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f $(OTHERDIR_MAIN)/*.mod
	@rm -f $(OTHERDIR_MOD)/*.mod
	@rm -f *.exe
	@rm -f *.mod

## clean_all: Same as "clean" and, additionally, the results and figures are also deleted (the repository returns to its original state).
.PHONY: clean_all
clean_all:
	@rm -f $(OTHERDIR_RESULTS)/*
	@rm -f $(OTHERDIR_FIG)/*
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f $(OTHERDIR_MAIN)/*.mod
	@rm -f $(OTHERDIR_MOD)/*.mod
	@rm -f *.dat
	@rm -f *.xyz
	@rm -f *.png
	@rm -f *.exe
	@rm -f *.mod


## The commands are explained in the order in which it is recommended to use them (if you want to use option 2).

.PHONY : help
help :
	@sed -n 's/^##//p' Makefile
