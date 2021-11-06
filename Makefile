all: sequential sequential_parallelo

sequential:
	mpicc -O3 -o sequential sequential.c
	
sequential_parallelo:
	mpicc -O3 -o sequential_parallelo sequential_parallelo.c
	
clean:
	rm -f sequential sequential_parallelo *.o