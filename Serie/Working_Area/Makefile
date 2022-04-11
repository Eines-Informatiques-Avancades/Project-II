SHELL := /bin/bash # necessary to use "time" shell command
# define variables
OTHERDIR_MAIN=../Code

OTHERDIR_MOD=../Code/Modules

OTHERDIR_RESULTS=../Results/Data

OTHERDIR_FIG=../Results/Figures

MOD_ORDENADOS=../Code/Modules/Constants_module.f90 ../Code/Modules/Pbc_module.f90 ../Code/Modules/Initial_state_module.f90 ../Code/Modules/Forces_module.f90 ../Code/Modules/Integrator_module.f90

OBJS= $(OTHERDIR_MAIN)/paralel.o $(OTHERDIR_MAIN)/constantes.o $(OTHERDIR_MAIN)/initial_state.o $(OTHERDIR_MAIN)/pbc.o $(OTHERDIR_MAIN)/forces.o $(OTHERDIR_MAIN)/integrator.o $(OTHERDIR_MAIN)/main.o

#----------------------------------------------------------------------------------------------------------------------------------#

## Information about how to use this makefile.
##
## Option 1: After entering the parameters in parameters.txt, you can use make all to realize the whole process: Simulation + Statistical analysis + Figures.
##           If you use this method, the data will automatically be saved in the folder Results/Data and the figures in the folder Results/Figures.
##           Additionally, the unnecessary files will be removed.
##
## Option 2: After entering the parameters in parameters.txt, you can use the commands that are explained below to realize the process step by step (less practical in my opinion).
##           If you want to do this, you must keep in mind that you only can analyse files that are in this folder, so do not use the move command until you have analysed the files.
##           Additionally, statistical analysis programs work with the files that are specified in parameters.txt, so make sure the file names are correct.
##
## Commands
##
## all       : compile fortran code, run the code, realize the statistical analysis and the figures, move the results to their folder and remove the useless files.
all: compile run statistics move clean

## compile   : compile the fortran code.
compile: paralel constantes initial_state pbc forces integrator main gofr program program_gofr

constantes: $(OTHERDIR_MOD)/Constants_module.f90
	@gfortran -c $(OTHERDIR_MOD)/Constants_module.f90 -o $(OTHERDIR_MAIN)/constantes.o

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
paralel: $(OTHERDIR_MOD)/paralel.f90
	@gfortran -c $(OTHERDIR_MOD)/paralel.f90 -o $(OTHERDIR_MAIN)/paralel.o 
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
initial_state: $(OTHERDIR_MOD)/Initial_state_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	@gfortran -c $(OTHERDIR_MOD)/Initial_state_module.f90 -o $(OTHERDIR_MAIN)/initial_state.o

pbc: $(OTHERDIR_MOD)/Pbc_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	@gfortran -c $(OTHERDIR_MOD)/Pbc_module.f90 -o $(OTHERDIR_MAIN)/pbc.o

forces: $(OTHERDIR_MOD)/Forces_module.f90 $(OTHERDIR_MOD)/Pbc_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	@gfortran -c $(OTHERDIR_MOD)/Forces_module.f90 -o $(OTHERDIR_MAIN)/forces.o	

integrator: $(OTHERDIR_MOD)/Integrator_module.f90 $(OTHERDIR_MOD)/Forces_module.f90 $(OTHERDIR_MOD)/Pbc_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	@gfortran -c $(OTHERDIR_MOD)/Integrator_module.f90 -o $(OTHERDIR_MAIN)/integrator.o

main: $(OTHERDIR_MAIN)/Main.f90 $(MOD_ORDENADOS)
	@mpif90 -g -o $(OTHERDIR_MAIN)/main.o $(OTHERDIR_MAIN)/Main.f90

gofr: $(OTHERDIR_MAIN)/Gofr.f90
	@gfortran -c $(OTHERDIR_MAIN)/Gofr.f90 -o $(OTHERDIR_MAIN)/gofr.o
	
program: $(OBJS)
	@mpirun -np 4 program.exe $(OBJS)

program_gofr: $(OTHERDIR_MAIN)/gofr.o
	@gfortran -o program_gofr.exe $(OTHERDIR_MAIN)/gofr.o

## run       : run the fortran code.
run: 
	@time ./program.exe
	@./program_gofr.exe

## statistics: realize the statistical analysis and the figures. You will be asked some questions.
.PHONY: statistics
statistics:
	@python $(OTHERDIR_MAIN)/generate_figures.py parameters.txt
	@python $(OTHERDIR_MAIN)/Block_avarage.py parameters.txt

## move      : move the results (files .dat and .png) to Results folder (Data and Figures).
.PHONY: move	
move: 
	@mv *.dat $(OTHERDIR_RESULTS)
	@mv *.png $(OTHERDIR_FIG)
	@mv *.xyz $(OTHERDIR_RESULTS)

## clean     : remove objects, executables and modules.
.PHONY: clean
clean:
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f *.exe
	@rm -f *.mod
	@rm -f gofr_params.txt

## clean_all : remove the same than clean and, additionally, the results and figures (the repository returns to its original state).
.PHONY: clean_all
clean_all:
	@rm -f $(OTHERDIR_RESULTS)/*
	@rm -f $(OTHERDIR_FIG)/*
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f gofr_params.txt
	@rm -f *.dat
	@rm -f *.xyz
	@rm -f *.png
	@rm -f *.exe
	@rm -f *.mod

# not available 
#.PHONY: vmd
#vmd:
#	vmd $(OTHERDIR_RESULTS)/vmd.xyz

##
## The commands are explained in the order in which it is recommended to use them (if you want to use option 2).

.PHONY : help
help :
	@sed -n 's/^##//p' Makefile
