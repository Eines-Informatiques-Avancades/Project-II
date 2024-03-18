<h1 align="center">
  <a name="logo" href="https://web.ub.edu/web/estudis/w/masteruniversitari-md308"><img src="https://web.ub.edu/documents/2710030/4734497/logo_UB_nou.jpg_ca.jpg/fbe1020f-8346-2e48-247d-cde35b5e38fc?t=1677492326543" alt="UB" width="200"></a>
  <br>
  Molecular Dynamics Simulation of a Van der Waals Gas
</h1>
<h4 align="center">Project-II - Eines Informàtiques Avançades</h4>
<p align="center"><a align="center" target="_blank" href="https://vcloudinfo.us12.list-manage.com/subscribe?u=45cab4343ffdbeb9667c28a26&id=e01847e94f"><img src="https://feeds.feedburner.com/RecentCommitsToBearStoneHA.1.gif" alt="Recent Commits to Bear Stone Smart Home" style="border:0"></a></p>
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

* <a href="https://github.com/Qbadosfe"><img src="https://avatars.githubusercontent.com/u/162143734?v=4" title="Qbadosfe" width="50" height="50"></a> **Quim Badosa |** Forces
* **Guillem Arasa |** Initial and boundary conditions
* **Rocío Aragoneses |**  Statistics and visualization
* **Emma Valdés |** I/O and integration
* **Paula Sierra |** Coordination

## Version History

* 0.1
    * Initial Release

