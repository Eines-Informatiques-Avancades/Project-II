# Compilator, el optimizator and flags variables
FC = mpif90
MC = mpirun
NP = 2
OPT = -O2
#FLAGS

# List of source files .f90
SRCS = initialization.f90 pbc.f90 verlet.f90 forces.f90 integrators.f90 readers.f90 main.f90

# List of object files .o
OBJS = $(SRCS:.f90=.o)

# Executable name
EXE = MD.exe

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
	mkdir -p results
	$(MC) -np $(NP) ./$(EXE)
	gnuplot visualization.gn
	mv *.dat *.png ./results

# Delete files .o, .mod and .out
clean:
	rm -r $(OBJS) *.mod *out *.err *.out HELLO

# Send process to queue iqtc04 or iqtc07
qsub4:
	@qsub -q iqtc04.q openmpi.sub
	@echo "Sent to queue iqtc04"

qsub7:
	@qsub -q iqtc07.q openmpi.sub
	@echo "Sent to queue iqtc07"

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
	@echo "  make qsub4  : Runs the process in the queue iqtc04"
	@echo "  make qsub7  : Runs the process in the queue iqtc07"
	@echo "  make clean  : Deletes files .o, .mod and .out"
	@echo "  make print  : Shows the Makefile variables used"
	@echo "  make help   : Shows this menu"
