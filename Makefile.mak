# Compilator, el optimizator and flags variables
FC = gfortran
OPT = -O3
FLAGS = -Wall -Wextra

# List of source files .f90
SRCS = $(wildcard *.f90)

# List of object files .o
OBJS = $(SRCS:.f90=.o)

# Executable name
EXE = calcular

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
clean:
	rm -f $(OBJS) *.mod

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
	@echo "  make -f Makefile.mak all    : Compiles the program $(EXE)"
	@echo "  make -f Makefile.mak run    : Runs the program $(EXE)"
	@echo "  make -f Makefile.mak clean  : Deletes files .o and .mod"
	@echo "  make -f Makefile.mak print  : Shows the Makefile variables used"
	@echo "  make -f Makefile.mak help   : Shows this menu"