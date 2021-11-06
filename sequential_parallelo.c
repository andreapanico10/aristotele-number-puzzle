#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


void permute(int *, int, int, int, int);
void swap(int*, int, int);
void printArray( int*, int);

int main(int argc, char **argv)
{
    int rank;
    int p;
    double elapsed_time;
    int i;
    int input;
    
    int *A;
      
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
   
    if(argc != 2){
      if (!rank)
           fprintf(stderr, "\n\nUsage: %s <max_number>\n\n", argv[0]);

      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    input = strtol(argv[1], NULL, 10);
    
    A = (int *) malloc(sizeof(int) * input);
    for(i = 1; i<=input;i++){
        A[i-1] = i;
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time =-MPI_Wtime();

    
    permute( A, 0, input, rank, p);

    
    elapsed_time += MPI_Wtime();
    if(!rank){
        fprintf(stdout,"Time: %lf\n",elapsed_time);
    }
    
    MPI_Finalize();
    
    free(A), A = NULL;

    return 0;


}

/* swap 2 values by index */
void swap(int* arr, int a, int b)
{
  int tmp = arr[a];
  arr[a] = arr[b];
  arr[b] = tmp;
}


/* print len elements */
void printArray( int* a, int len)
{
 size_t i = 0;
  
 fprintf(stdout, "[ " );

 for ( i=0; i< len; i++) fprintf(stdout, "%d " , a[i] );

 fprintf(stdout, "]\n" );
}



/* permute an array recursively */
void permute(int *arr, int start, int end, int rank, int p)
{
    if(arr[0] <= (rank+1)*end/p){
        int i;

        if(start == end) /* this function is done */
        {
         //printArray(arr, end);
         return;
        }

        permute(arr, start + 1, end, rank, p); /* start at next element */

        /* permute remaining elements recursively */
        for(i = start + 1; i < end; i++)
        {
           if( arr[start] == arr[i] ) continue; /* skip */

           swap(arr, start, i);

           permute(arr, start + 1, end, rank, p);

           swap(arr, start, i); /* restore element order */

        }

    }
}
 


