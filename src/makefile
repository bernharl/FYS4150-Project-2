CPPflags= c++ -O3

LIB = -larmadillo

TESTPROG= testcode.out
MAINPROG= mainprog.out

all:	${TESTPROG} ${MAINPROG}

${TESTPROG} :	   		tests-main.o  test-functions.o jacobi_algorithm.o
					${CPPflags} tests-main.o test-functions.o jacobi_algorithm.o ${LIB} -o ${TESTPROG}

tests-main.o :			tests-main.cpp 
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

jacobi_algorithm.o :		jacobi_algorithm.cpp 
		        	${CPPflags} -c jacobi_algorithm.cpp


${MAINPROG}:	main.cpp
				${CPPflags} main.cpp ${LIB} -o ${MAINPROG}



