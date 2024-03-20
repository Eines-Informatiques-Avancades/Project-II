# Compilator, el optimizator and flags variables
FC = gfortran
OPT = -O3
FLAGS = -Wall -Wextra

# List of source files .f90
SRCS = inicialitzacio.f90 pbc_EIA.f90 potential.f90 integrators.f90 readers_mod.f90 main.f90

# List of object files .o
OBJS = $(SRCS:.f90=.o)

# Executable name
EXE = a.out

# Files compilation
all: $(EXE)

# Executable file compilation from the object files
$(EXE): $(OBJS)
	$(FC) $(OPT) $(FLAGS) -o $@ $^

# Source file compilation from .f90 to .o
%.o: %.f90
	$(FC) $(OPT) $(FLAGS) -c $<

# Run 
run: $(EXE)
	./$(EXE)

# Delete files .o and .mod
tidy:
	rm -f $(OBJS) *.mod

# Delete files .o and .mod
clean:
	rm -f $(OBJS) *.mod *.dat *.out

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
	@echo "  make tidy   : Deletes files .o and .mod"
	@echo "  make clean  : Deletes files .o, .mod, .dat and .out"
	@echo "  make print  : Shows the Makefile variables used"
	@echo "  make help   : Shows this menu"