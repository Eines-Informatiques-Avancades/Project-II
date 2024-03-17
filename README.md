# Project-II

## Introduction

This project consists on a simple molecular dynamics program to simulate a Van der Waals gas of particles. It allows to obtain results indicating the evolution of the energy (kinetic, potential, total), temperature and pressure.

## Description

The different modules that can be found in this Repo are:

* Initialization: Method that will set all the initial particles within a periodical box and assign their initial state.
* Boundary conditions: Method to assign particle coordinates within periodical box.
* Forces: Method that calculates the van der Waals forces at each particle.
* Integrators: Method to integrate Newton’s equations.
* Statistics: Method to determine all averaged values of selected variables during simulation trajectory and showing its final values at the end of the trajectory.
* Visualization of results: post-processing method that will deal with obtained data from a given trajectory and will allow to show evolution of the different variables involved in the molecular dynamics and other statistics at the highest level.

## Getting Started

### Dependencies

- Gnuplot
- Fortran

### Installing

1. Download all the files with extension `.f90` and `.mak` and place them in the same directory.
2. Modify the compilator, optimizator and flag variables in `Makefile.mak` to your preference.
3. Modify the parameters of your simulation in `inputs.txt`.

### Executing program

Run the `Makefile.mak` script.
```
make -f Makefile.mak
```

## Help

The `Makefile.mak` script has a help menu that can be accessed running the following command:
```
make help
```

## Authors

* **Quim Badosa |** Forces
* **Guillem Arasa |** Initial and boundary conditions)
* **Rocío Aragoneses |**  Statistics and visualization
* **Emma Valdés |** I/O and integration
* **Paula Sierra |** Coordination

## Version History

* 0.1
    * Initial Release

