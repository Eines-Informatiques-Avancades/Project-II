<h1 align="center">
  <a name="logo" href="https://web.ub.edu/web/estudis/w/masteruniversitari-md308"><img src="https://web.ub.edu/documents/2710030/4734497/logo_UB_nou.jpg_ca.jpg/fbe1020f-8346-2e48-247d-cde35b5e38fc?t=1677492326543" alt="UB" width="200"></a>
  <br>
  Molecular Dynamics Simulation of a Van der Waals Gas
</h1>
<h4 align="center">Project-II - Eines Informàtiques Avançades</h4>
<div align="center">

  <h4>
    <a href="https://github.com/Eines-Informatiques-Avancades/Project-II/commits/master"><img src="https://img.shields.io/github/last-commit/Eines-Informatiques-Avancades/Project-II.svg?style=plasticr"/></a>
        <a href="https://github.com/Eines-Informatiques-Avancades/Project-II/commits/master"><img src="https://img.shields.io/github/commit-activity/y/Eines-Informatiques-Avancades/Project-II.svg?style=plasticr"/></a>

  </h4>
</div>
<p><font size="3">

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

1. Run the `Makefile.mak` script.
```
make -f Makefile.mak
```
2. Execute the output program.
```
./calcular
```

## Help

The `Makefile.mak` script has a help menu that can be accessed running the following command:
```
make -f Makefile.mak help
```

## Authors

* <a href="https://github.com/Qbadosfe"><img src="https://avatars.githubusercontent.com/u/162143734?v=4" title="Qbadosfe" width="25" height="25"></a> **Quim Badosa |** Forces
* <a href="https://github.com/guillemares"><img src="https://avatars.githubusercontent.com/u/144935605?v=4" title="guillemares" width="25" height="25"></a> **Guillem Arasa |** Initial and boundary conditions
* <a href="https://github.com/Qbadosfe"><img src="https://avatars.githubusercontent.com/u/162143734?v=4" title="Rocio" width="25" height="25"></a> **Rocío Aragoneses |**  Statistics and visualization
* <a href="https://github.com/evaldesmartin"><img src="https://avatars.githubusercontent.com/u/125901837?v=4" title="evaldesmartin" width="25" height="25"></a> **Emma Valdés |** I/O and integration
* <a href="https://github.com/psierrva"><img src="https://avatars.githubusercontent.com/u/162144063?v=4" title="psierrva" width="25" height="25"></a> **Paula Sierra |** Coordination

## Version History

* 0.1
    * Initial Release

