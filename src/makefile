# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we need
LIB = -larmadillo
# Here we define the name of the executable
PROG= testcode
${PROG} :	   	tests-main.o  test-functions.o main.o
			${CPPflags} tests-main.o test-functions.o main.o ${LIB} -o ${PROG}

main.o :		main.cpp 
		        	${CPPflags} -c main.cpp

tests-main.o :			tests-main.cpp 
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp
