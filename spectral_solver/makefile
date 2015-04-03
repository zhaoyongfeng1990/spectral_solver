solver: derivative.o fun.o initialization.o iteration.o main.o odesolver.o printstatus.o solver.o
	g++-4.9 *.o -o solver -O3 -lm -lgsl -lgslcblas -lfftw3 -lfftw3_threads -fopenmp
derivative.o: derivative.cpp solver.h parameters.h
	g++-4.9 -O3 -c derivative.cpp -std=gnu++11
solver.o: solver.cpp solver.h parameters.h
	g++-4.9 -O3 -c solver.cpp -std=gnu++11
fun.o: fun.cpp solver.h parameters.h
	g++-4.9 -O3 -c fun.cpp -std=gnu++11
odesolver.o: odesolver.cpp solver.h parameters.h
	g++-4.9 -O3 -c odesolver.cpp -std=gnu++11
iteration.o: iteration.cpp solver.h parameters.h
	g++-4.9 -O3 -c iteration.cpp -std=gnu++11
initialization.o: initialization.cpp solver.h parameters.h
	g++-4.9 -O3 -c initialization.cpp -std=gnu++11
printstatus.o: printstatus.cpp solver.h parameters.h
	g++-4.9 -O3 -c printstatus.cpp -std=gnu++11
main.o: main.cpp solver.h parameters.h
	g++-4.9 -O3 -c main.cpp -std=gnu++11
test:
	make solver
	./solver
clean:
	rm ./*.o
remove:
	rm ./solver
	rm ./*.txt
