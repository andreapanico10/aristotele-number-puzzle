#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    int rank;
    int p;
    double elapsed_time;
    int i;
        
    int A[] = {10,2,10,10,10,10,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,10,10,10,10,10,10,10,10,10,10,10,5,3,4,23,4,1,1,1,1,1,1,1,1,1};
      
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time =-MPI_Wtime();

    
    
    
    elapsed_time += MPI_Wtime();

    
    MPI_Finalize();

    return 0;

}


