all: sequential 

sequential:
	mpicc -O3 -o sequential sequential.c
	
clean:
	rm -f sequential *.o