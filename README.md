# FYS4150-Project-2
Second project in Computational Physics(FYS4150)

## Build Instructions:

* Run main_script.sh, this will build the source code, run tests and the main executable, generate plots using the Python library Matplotlib, and (optional) build the latex report for the project.
* Should you want to run compile the c++ source code manually, run "make" in src/, this should generate two executables, one for unit tests using Catch2 (https://github.com/catchorg/Catch2) and one for running the main calculations.
&nbsp;

## Run Instructions:
* src/mainprog.out Is compiled when running make, this executable takes two input arguments: n: dimension of nxn matrix, eps: tolerance for treating off-diagonal elements as zero in the Jacobi Algorithm, ![equation](http://www.sciweavers.org/tex2img.php?eq=%24%5Crho_n%24&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0): What is considered "infinity" in the quantum system.
* src/testcode.out Run this for unit tests
&nbsp;

## Structure: 

* src/Jacobi_argorithm.cpp Contains all methods used in calculations.
* src/main.cpp Runs the methods in Jacobi_algorithm.cpp. Results are saved in src/eigendata.dat
* src/tests_main.cpp Contains all unit tests run by Catch2.
* src/plot_infinity.py Plots the absolute errors of the first four eigenvalues of quantum dots in three dimensions with a single electron. Plots are saved in /doc/figures
* src/plot_states.py Plots the ground state eigenvectors for a harmonic oscillator well with repulsive Coulomb interaction. Plots are saved in doc/figures/
&nbsp;
 
