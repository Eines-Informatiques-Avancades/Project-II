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

This project consists on a simple molecular dynamics program to simulate a Van der Waals gas of particles. It allows to obtain results indicating the evolution of the energy (kinetic, potential, total), temperature, pressure and radial distribution function.


## Description

The different modules that can be found in this Repo are:

* **Initialization:** Method that will set all the initial particles within a periodical box and assign their initial positions and velocities.
* **Boundary conditions:** Method to assign particle coordinates within periodical box.
* **Forces:** Method that calculates the van der Waals forces at each particle.
* **Integrators:** Method to integrate Newton’s equations. 
* **Statistics:** Method to determine all averaged values of selected variables during simulation trajectory and showing its final values at the end of the trajectory.
* **Visualization of results:** post-processing method that will deal with obtained data from a given trajectory and will allow to show evolution of the energy, temperature, pressure and radial distribution of the system.

## Getting Started

### Dependencies

- Fortran (Recommended gfortran)
- OpenMPI (Version 3.1.3_ics-2015.0 or newer)
- Gnuplot

### Installing and Preparation

1. Clone this repository or download raw `.zip` file in your local directory. For both Serial or Parallel execution, make sure the files remain in the original directory.

2. Modify the compilator, optimizator and flag variables in `Makefile` to your preference.

   **_Note:_** _The default compilator and optimizator are `gfortran` and `-O3` for **Serial**._

   **_Note:_** _The default compilator, executer and optimizator are `mpif90`, `mpirun` and `-O2` for **Parallel**. The default number of processors is `4`_.
   
3. Modify the parameters of your simulation in `parameters.nml`.

   **_Note:_** _The default parameters relate to a system of the noble gas Kr at 300K_ (Rutkai, G. et al, 2016. doi: 10.1080/00268976.2016.1246760).

### Executing program in Serial:

1. Run the `Makefile` script.
```
make -f Makefile
```
2. Execute the output program by running:
```
make run
```
This also generates the final output plots and data in a `results/` folder. Alternatively, you can run the program with `./a.out`, and generate your plots writing `gnuplot visualization.gn`.

### Executing program in Parallel:

1. Run the `Makefile` script.
```
make -f Makefile
```
2. **For local execution:** Execute the output program by running:
```
make run
```
This also generates the final output plots and data in a `results/` folder. Alternatively, you can run the program with `MD.exe`, and generate your plots writing `gnuplot visualization.gn`.

3. **For cluster execution:** Execute the output program by running:
```
make qsub4 # for execution in iqtc04.q
make qsub7 # for execution in iqtc07.q
```

> Before the execution, make sure you change the environment of your cluster in `openmpi.sub` by choosing the modules required.


This also generates the final output plots and data in a `HELLO/` folder. A `hello.log`, `hello.out` and `hello.err`. Alternatively, you can run the program with `MD.exe`, and generate your plots writing `gnuplot visualization.gn`.

## Help

The `Makefile` script has a help menu that can be accessed running the following command:
```
make help
```
An example of the results obtained for a system of Kr gas at 300K are presented below

<a name="Energies" href="https://github.com/Eines-Informatiques-Avancades/Project-II/tree/master/Results"><img src="https://github.com/Eines-Informatiques-Avancades/Project-II/blob/master/Results/ResultsParallel/energy_plot.png?raw=true" alt="Energies" width="400"></a>
<a name="Temperature" href="https://github.com/Eines-Informatiques-Avancades/Project-II/tree/master/Results"><img src="https://github.com/Eines-Informatiques-Avancades/Project-II/blob/master/Results/ResultsParallel/temp_plot.png?raw=true" alt="Temperature" width="400"></a>
<a name="Pressure" href="https://github.com/Eines-Informatiques-Avancades/Project-II/tree/master/Results"><img src="https://github.com/Eines-Informatiques-Avancades/Project-II/blob/master/Results/ResultsParallel/pressure_plot.png?raw=true" alt="UPC" width="400"></a>
<a name="g(r)" href="https://github.com/Eines-Informatiques-Avancades/Project-II/tree/master/Results"><img src="https://github.com/Eines-Informatiques-Avancades/Project-II/blob/master/Results/ResultsParallel/g_r.png?raw=true" alt="g(r)" width="400"></a>

In addition, the folder `Results` contains more information about this example. Remember that you can adjust the parameters to your simulation modifying the file `parameters.nml`.

* Values of Energy, Pressure, Instant Temperature, Radial Distribution in `.dat` format.
* Plots of Energy, Pressure, Instant Temperature, Radial Distribution in `.png` format.
* Trajectory file of the particles for VMD analysis in `.xyz` format.

  
## Authors

* <a href="https://github.com/Qbadosfe"><img src="https://avatars.githubusercontent.com/u/162143734?v=4" title="Qbadosfe" width="25" height="25"></a> **Quim Badosa |** Forces
* <a href="https://github.com/guillemares"><img src="https://avatars.githubusercontent.com/u/144935605?v=4" title="guillemares" width="25" height="25"></a> **Guillem Arasa |** Initial and boundary conditions
* <a href="https://github.com/rocii389"><img src="https://avatars.githubusercontent.com/u/150451845?v=4" title="rocii389" width="25" height="25"></a> **Rocío Aragoneses |**  Statistics and visualization
* <a href="https://github.com/evaldesmartin"><img src="https://avatars.githubusercontent.com/u/125901837?v=4" title="evaldesmartin" width="25" height="25"></a> **Emma Valdés |** I/O and integration
* <a href="https://github.com/psierrva"><img src="https://avatars.githubusercontent.com/u/162144063?v=4" title="psierrva" width="25" height="25"></a> **Paula Sierra |** Coordination

## Version History

* 0.1
    * Initial Release

