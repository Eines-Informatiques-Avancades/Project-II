# Compilator, el optimizator and flags variables
FC = $(HOME)/opt/usr/local/bin/mpif90
MC = $(HOME)/opt/usr/local/bin/mpirun -np 4

# List of source files .f90
SRCS = initialization.f90 pbc.f90 forces.f90 integrators.f90 readers.f90 main.f90

# List of object files .o
OBJS = $(SRCS:.f90=)

# Executable name
EXE = a.out

run:
	$(FC) -g $(SRCS)
#	$(MC)

# Run 
#run: $(EXE)
#	mkdir -p results
#	./$(EXE)
#	gnuplot visualization.gn
#	mv *.dat *.png ./results

# Delete files .o, .mod and .out
clean:
	rm -f $(OBJS) *.mod *out

# See the variables used
print:
	@echo "FC = $(FC)"
	@echo "OPT = $(OPT)"
	@echo "FLAGS = $(FLAGS)"
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "EXE = $(EXE)"

# Help menu of the Makefile
help:
	@echo "Available options:"
	@echo "  make all    : Compiles the program $(EXE)"
	@echo "  make run    : Runs the program $(EXE)"
	@echo "  make clean  : Deletes files .o, .mod and .out"
	@echo "  make print  : Shows the Makefile variables used"
	@echo "  make help   : Shows this menu"