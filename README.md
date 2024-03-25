<h1 align="center">
  <a name="UB" href="https://web.ub.edu/web/estudis/w/masteruniversitari-md308"><img src="https://web.ub.edu/documents/2710030/4734497/logo_UB_nou.jpg_ca.jpg/fbe1020f-8346-2e48-247d-cde35b5e38fc?t=1677492326543" alt="UB" width="160"></a>
  <a name="UPC" href="https://web.ub.edu/web/estudis/w/masteruniversitari-md308"><img src="https://www.rediprogram.eu/wp-content/uploads/2021/08/politecnica-catalunya-300x300.png" alt="UPC" width="160"></a>
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

<p align="center">
  <a href="#introduction">Introduction</a> •
  <a href="#description">Description</a> •
  <a href="#getting-started">Getting Started</a> •
  <a href="#help">Help</a> •
  <a href="#authors">Authors</a> •
  <a href="#version-history">Version History</a>
</p>

## Introduction

This project consists on a simple molecular dynamics program to simulate a Van der Waals gas of particles. It allows to obtain results indicating the evolution of the energy (kinetic, potential, total), temperature and pressure.

## Description

The different modules that can be found in this Repo are:

* **Initialization:** Method that will set all the initial particles within a periodical box and assign their initial positions and velocities.
* **Boundary conditions:** Method to assign particle coordinates within periodical box.
* **Forces:** Method that calculates the van der Waals forces at each particle.
* **Integrators:** Method to integrate Newton’s equations. 
* **Statistics:** Method to determine all averaged values of selected variables during simulation trajectory and showing its final values at the end of the trajectory.
* **Visualization of results:** post-processing method that will deal with obtained data from a given trajectory and will allow to show evolution of the energy, temperature and pressure of the system.

## Getting Started

### Dependencies

- Fortran
- Gnuplot

### Installing

1. Download the folder Serial code. During all the execution, make sure the files with extension `.f90`, `.gn` and the `Makefile` remain in the same directory.
2. Modify the compilator, optimizator and flag variables in `Makefile` to your preference.

   **_Note:_** _The default compilator and optimizator are `gfortran` and `-O3`._
   
3. Modify the parameters of your simulation in `parameters.nml`.

   **_Note:_** _The default parameters relate to a system of Kr gas at 300K._

### Executing program

1. Run the `Makefile` script.
```
make -f Makefile
```
2. Execute the output program by running
```
make run
```
This also generates the final output plots. Alternatively, you can run the program with `./a.out`, and generate your plots writing `gnuplot visualization.gn`.

## Help

* The `Makefile` script has a help menu that can be accessed running the following command:
```
make help
```
* The folder `Results` contains an example of the results obtained for a system of Kr gas at 300K. Remember that you can adjust the parameters to your simulation modifying the file `parameters.nml`.
  
## Authors

* <a href="https://github.com/Qbadosfe"><img src="https://avatars.githubusercontent.com/u/162143734?v=4" title="Qbadosfe" width="25" height="25"></a> **Quim Badosa |** Forces
* <a href="https://github.com/guillemares"><img src="https://avatars.githubusercontent.com/u/144935605?v=4" title="guillemares" width="25" height="25"></a> **Guillem Arasa |** Initial and boundary conditions
* <a href="https://github.com/rocii389"><img src="https://avatars.githubusercontent.com/u/150451845?v=4" title="rocii389" width="25" height="25"></a> **Rocío Aragoneses |**  Statistics and visualization
* <a href="https://github.com/evaldesmartin"><img src="https://avatars.githubusercontent.com/u/125901837?v=4" title="evaldesmartin" width="25" height="25"></a> **Emma Valdés |** I/O and integration
* <a href="https://github.com/psierrva"><img src="https://avatars.githubusercontent.com/u/162144063?v=4" title="psierrva" width="25" height="25"></a> **Paula Sierra |** Coordination

## Version History

* 0.1
    * Initial Release

