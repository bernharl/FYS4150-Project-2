# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= testcode
${PROG} :	   	tests-main.o  test-functions.o jacobi_algorithm.o
			${CPPflags} tests-main.o test-functions.o jacobi_algorithm.o ${LIB} -o ${PROG}

tests-main.o :			tests-main.cpp 
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

jacobi_algorithm.o :		jacobi_algorithm.cpp 
		        	${CPPflags} -c jacobi_algorithm.cpp
