# FYS4150-Project-2
Second project in Computational Physics(FYS4150)

## Build Instructions:

* Run main_script.sh, this will build the source code, run tests and the main executable, generate plots using the Python library Matplotlib, and (optional) build the latex report for the project.
* Should you want to run compile the c++ source code manually, run "make" in src/, this should generate two executables, one for unit tests using Catch2 (https://github.com/catchorg/Catch2) and one for running the main calculations.
&nbsp;

## Structure: 

* src/Jacobi_argorithm.cpp Contains all methods used in calculations.
* src/main.cpp Runs the methods in Jacobi_algorithm.cpp.
* src/tests_main.cpp Contains all unit tests run by Catch2.
* src/plot_states.py Plots the ground state eigenvectors for a harmonic oscillator well with repulsive Coulomb interaction.
&nbsp;
 
