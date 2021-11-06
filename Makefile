all: sequential sequential_parallelo solver_github

sequential:
	mpicc -O3 -o sequential sequential.c
	
sequential_parallelo:
	mpicc -O3 -o sequential_parallelo sequential_parallelo.c

solver_github:
	mpicc -O3 -o solver_github solver_github.c
	
clean:
	rm -f sequential sequential_parallelo solver_github *.o