solver: derivative.cpp fun.cpp initialization.cpp iteration.cpp list.txt main.cpp odesolver.cpp parameters.h printstatus.cpp solver.cpp solver.h
	g++-4.9 -O3 -lm -lgsl -lgslcblas -lfftw3 -lfftw3_threads -fopenmp *.cpp -o solver

mpi: derivative.cpp fun.cpp initialization.cpp iteration.cpp list.txt main.cpp odesolver.cpp parameters.h printstatus.cpp solver.cpp solver.h
	mpic++ -O3 -lm -lgsl -lgslcblas -lfftw3 -lfftw3_threads *.cpp -o solver

run:
	mpirun -np 8 ./solver
test:
	make mpi
	make run

clean:
	rm ./solver
	rm ./*.txt
