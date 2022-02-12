CC=mpicc
SRC=.c MyMPI.c
EXEC=magic_hexagon
CFLAGS=-O3

all: $(EXEC)

magic_hexagon:
	$(CC) $(CFLAGS) -o $@ magic_hexagon.c 

clean:
	rm -f $(EXEC) *.o *~