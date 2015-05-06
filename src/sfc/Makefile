
test: evalhilbert.o hilbert.o testF.o
	gfortran -o test testF.o evalhilbert.o hilbert.o 

testF.o: testF.f95
	gfortran -c testF.f95

evalhilbert.o: evalhilbert.c
	gcc -c evalhilbert.c -I./

hilbert.o : hilbert.c
	gcc -c hilbert.c -I./
