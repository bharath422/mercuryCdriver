mercury_c: mercury.o main.o
	gfortran -o mercury_c mercury.o main.o

mercury_f: mercury.o main_mercury.o
	gfortran -o mercury_f mercury.o main_mercury.o

main_mercury.o: main_mercury.f
	gfortran -c main_mercury.f

mercury.o: mercury.f
	gfortran -c mercury.f

main.o: main.c
	gcc -c main.c

clean:
	rm -f mercury_c mercury_f *.o

cout:
	rm -f *.tmp *.dmp *.out

cleanall:
	make clean
	make cout
