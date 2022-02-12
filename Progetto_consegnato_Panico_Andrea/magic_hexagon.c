#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 15   /*row of input matrix: 15 constraints*/
#define M 20   /*columns of input matrix: 19 variables + 1 column of known terms*/

typedef enum { false, true } bool;
double elapsed_time;  /* Time to find, count solutions */
int n_solutions = 0;  /* number of solutions founded (needed only due to print purpouses)*/
int founded = 0;      /* it is a value 0/1. When it is 1 for 1 process, the search will end*/
int founder = 0;      /* at the end of execution all the processes know the process-id of the founder*/
int *empty_diagonal;  /* array that own the indexes of the blank row in the system in gaussian form 19x19; we can forecast the len of *empty_diagonal (=7 = 19 variables - 12 independent equations)*/
int id;               /* Process rank */
int p;                /* Number of processes */
int iterations = 0;   /* count the number of iterations in the recursive computation of all the permutations*/
int analyzed ;        /* is the number of the last iteration analyzed by a process*/
int gaussian[M-1][M]; /* matrix 19*20   |A||B| B is juxtaposed to A */
int A[N][M];          /* input matrix 15*20 that contains the problem constraints*/
int **solutions;

bool back_substitution(int matrix[M-1][M]);
void get_final_values (int gaussian[M-1][M]);
void find_all_solutions(int *array);
bool validation(int *array);
void swap_int(int *a, int *b);
void permutation(int *a, int l, int r);
void add_solution(int *solution);
void valid_permutation(int *permutation, int size);
void transform_array_in_rotation(int *array, int *rotation);
void transform_rotation_in_array(int *array, int *rotation);
void find_flip(int *array);
void get_combination(int arr[], int n,int r);
void combination_util(int arr[], int data[], int start, int end,int index, int r);
int read_matrix(int rows, int cols, int (*a)[cols], const char* filename);
void print_hexagons();
void swap_row(int idx1, int idx2);
void multiply(int i, int c);
void add(int i, int j, int c);

int main (int argc, char *argv[]) {

    /* Initialize MPI*/
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    /*domain has domain values*/
    int domain[M-1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    int r = 7 ; /*number of unknown variables = 19 variables - 12 independent equations*/

    /*Getting matrix A of input constraints*/
    read_matrix(N,M, A,"initial_matrix.txt" );
    /*  a b c d e f g h i j k l m n o p q r s
       -----------------------------------------
        1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 38
        0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 38
        0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 38
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 38
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 38
        1 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 38
        0 1 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 38
        0 0 1 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0 0 38
        0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0 38
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 38
        0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 38
        0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 38
        1 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 38
        0 1 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 38
        0 0 1 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 38
      */

    int row_index = 0;
    int col_index = 0;
    int broken = 0;

    /* Start timer */
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();

    /*the following while loop transform the input matrix 15*20 in a matrix 15*20 where the 12 equations are independent + 3 blank rows (gaussian elimination) */
    while (row_index < N && col_index < M-1){
        int j;
        for ( j = row_index; j < N; ++j) {
            broken = 0; /* keep track if the for loop is broken*/
            if(A[j][col_index]){
                swap_row(row_index,j);
                broken = 1;
                break;
            }
        }
        /*if the for loop is not broken, increment the column_index*/
        if ( broken == 0 ){
            col_index = col_index+ 1;
            continue;
        }

        /* row is a pointer to the first element of row_index A's row */
        int *row;
        row = A[row_index];
        int leading_value = row[col_index];
        multiply(row_index, 1 / leading_value);

        /*for loop update the matrix adding the changed row -> row_index to all the other rows*/
        for (int j = 0; j < N; ++j) {
                if (j == row_index)
                    continue;
                else
                    add(j, row_index, -A[j][col_index]);
        }

        row_index = row_index + 1;
        col_index = col_index + 1;
    }

    int k = 0;
    int displacement = 0; /*meaning explained in the following rows*/

    /*how to generate the final gaussian elimination matrix in 19*20 with 7 blank rows (now we are in a situation 12*20 with 3 blank row , must be added the 4 blank rows and put the total 7 blank rows in the right places)*/
    /*
    1 0 0 0 0 0 0 0 0  1  1 0 0  1  2  1 0  1  1  76                    1 0 0 0 0 0 0 0 0  1  1 0 0  1  2  1 0  1  1  76
    0 1 0 0 0 0 0 0 0 -1  0 0 0 -1 -1  0 0  0  0  0                     0 1 0 0 0 0 0 0 0 -1  0 0 0 -1 -1  0 0  0  0  0
    0 0 1 0 0 0 0 0 0  0 -1 0 0  0 -1 -1 0 -1 -1 -38                    0 0 1 0 0 0 0 0 0  0 -1 0 0  0 -1 -1 0 -1 -1 -38
    0 0 0 1 0 0 0 0 0 -1 -1 0 0  0 -1  0 0  0  0  0                     0 0 0 1 0 0 0 0 0 -1 -1 0 0  0 -1  0 0  0  0  0
    0 0 0 0 1 0 0 0 0  0 -1 0 0 -1 -1 -1 0 -1  0 -38                    0 0 0 0 1 0 0 0 0  0 -1 0 0 -1 -1 -1 0 -1  0 -38
    0 0 0 0 0 1 0 0 0  1  1 0 0  1  1  1 0  0  0  38                    0 0 0 0 0 1 0 0 0  1  1 0 0  1  1  1 0  0  0  38
    0 0 0 0 0 0 1 0 0  0  1 0 0  0  1  0 0  1  0  38                    0 0 0 0 0 0 1 0 0  0  1 0 0  0  1  0 0  1  0  38
    0 0 0 0 0 0 0 1 0  0  0 0 0 -1 -1 -1 0 -1 -1 -38                    0 0 0 0 0 0 0 1 0  0  0 0 0 -1 -1 -1 0 -1 -1 -38
    0 0 0 0 0 0 0 0 1  1  1 0 0  1  1  0 0  1  0  38                    0 0 0 0 0 0 0 0 1  1  1 0 0  1  1  0 0  1  0  38
    0 0 0 0 0 0 0 0 0  0  0 1 0  0  0  1 0  0  1  38    ---> to -->     0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0  0
    0 0 0 0 0 0 0 0 0  0  0 0 1  1  1  1 0  0  0  38                    0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0  0
    0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 1  1  1   38                   0 0 0 0 0 0 0 0 0  0  0 1 0  0  0  1 0  0  1  38
    0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0                    0 0 0 0 0 0 0 0 0  0  0 0 1  1  1  1 0  0  0  38
    0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0                    0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0
    0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0                    0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0
                                                                        0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0
                                                                        0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 1  1  1   38
                                                                        0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0
                                                                        0 0 0 0 0 0 0 0 0  0  0 0 0  0  0  0 0  0  0   0
     */
    /*for loop look ideally at the element on the main diagonal, but due to the not squared matrix, when the diagonal element is 0
     * it must be added a blank row; example:
     *  0 0 0 0 0 0 0 0 0  0  0 1 0  0  0  1 0  0  1  38 is the 9-th row of A matrix (9 starting from index 0)
     * The 9-th element ☝🏻           is 0 -> if violation
     * The 1st element !0 is    ☝     (2 position after 9 -> so 2 blank rows must be added )

     *
     */


    int elements = 0; //keep track of the index of the variables unassigned in the gaussian form;
    empty_diagonal= (int*) malloc((elements+1) * sizeof(int));

    for (int i = 0; i < M-1; ++i) {
        if (A[k][k + displacement]){
            for (int j = 0; j < M; ++j) {
                gaussian[i][j] = A[k][j];
            }
            k ++;
        }
        else {
            displacement++;
            /* adding i index in the list of blank row in 19x20 the gaussian matrix */
            elements++;
            empty_diagonal = realloc(empty_diagonal, (elements+1) * sizeof(int));
            empty_diagonal[elements - 1 ] = i;
        }
    }
    /*any process has a different value for analyzed, we start from the value -p + id in order to obtain round robin work assignment*/
    analyzed = -p +id;

    /*Main work
     * LOGICAL FLOW
     *  1) get_combination call combination_util
     *  2) combination_util generate a combination and call permutation
     *  3) permutation produces all the permutations of the combination and call valid_permutation and call back_substitution
     *  4) back_substitution take the 7 elements of the permutation and complete the 12x20 gaussian matrix to a 19*20 matrix, validate the domain constraints and call
     *  5) if back_substitution is correct, get_final_values is called that extract the 19 known terms of the variables and it calls find_all_solutions
     *  6) find_all_solutions find the solutions linked to the found one (rotations and flippings)
     *  7) communicate to the other processes that the solutions were found; come back to main and print_hexagons()
     */


    get_combination(domain, M-1, r);

    /*End time*/
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    /*All the process know founder id, but only the founder has got the 12 solutions,so he print them*/
    if(id == founder){
        print_hexagons();
        fprintf (stdout,"Execution time %8.6f\n\n", elapsed_time);
        fprintf(stdout,"\nProcess %d found the solutions\n",id);
        FILE *out_file = fopen("times.txt", "a"); // append
        if(!out_file){
            fprintf(stderr,"Error in writing output in file\n");
        }
        fprintf(out_file, "%8.6f\n", elapsed_time); // write to file
        fclose(out_file);


    } else
        fprintf(stdout,"Proc %d  leaves \n",id);

    MPI_Finalize();
     free(solutions), solutions = NULL;
     free(empty_diagonal), empty_diagonal = NULL;
    return 0;
}

/* function that get all combinations of size r (in our case r = 7 -> the seven unknown variables)
  in input_array[] of size n (in our case n = 19). This function mainly uses combinationUtil() */
void get_combination(int input_array[], int n, int r)
{
    // A temporary array to store all combination one by one
    int data[r];
    // Print all combination using temporary array 'data[]'
    combination_util(input_array, data, 0, n-1, 0, r);
}

/* arr[]  ---> Input Array
   data[] ---> Temporary array to store current combination
   start & end ---> Starting and Ending indexes in arr[]
   index  ---> Current index in data[]
   r ---> Size of a combination to be printed */
void combination_util(int arr[], int data[], int start, int end,int index, int r){
    // Current combination is ready to be printed, print it

        if (index == r) {
            if (!founded) {
                permutation(data, 0, 6);
            } else {
                return;
            }
        }
        /* replace index with all possible elements. The condition
         * "end-i+1 >= r-index" makes sure that including one element
         * at index will make a combination with remaining elements
         * at remaining positions
         */

        for (int i = start; i <= end && end - i + 1 >= r - index; i++) {
            data[index] = arr[i];
            combination_util(arr, data, i + 1, end, index + 1, r);
        }
}


// this method find upper triangular matrix' n solutions
bool back_substitution(int matrix[M-1][M]){

    for (int i = (M-1)-1; i >= 0 ; --i) {
        double value = matrix[i][M-1]/matrix[i][i];
        //excluding all the solutions xi that have value b <= 0 || > 19 [constraint violation]
        if(value >19.0 || value < 1.0)
            return false;

        for (int j = 0; j < i; ++j) {
            if(matrix[j][i] != 0){
                matrix[j][M-1] = matrix[j][M-1] - (value * matrix[j][i]);
                matrix[j][i] = 0;
            }
        }
    }
    return true;
}

void get_final_values(int gaussian[M-1][M]){

    int *array =  (int*) malloc((M-1) * sizeof(int));
    if(!array){
        fprintf(stderr,"Malloc error\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    int index = 0;

    // array is the list of the known terms,i.e. the last column
    for (int i = 0; i < M; ++i) {
        array[index] = gaussian[i][M-1];
        index++;
    }

    //calling the validation function that check the constraints
    bool valid = validation(array);

    if(valid == true){
        //if there is no violation -> find all the solutions starting from the first founded and comunicate to the other processses
        founded = 1;
        founder = id;
        find_all_solutions(array);

        int next = (id+1) % p;
        // but we fix the problematic case -> next = 0
        if (id == p-1)
            next = 0;
        int message = founder;
        /*Send to the successor in the ring network the founder-id (In this way the receiver known that a solution was found and who found it)*/
        MPI_Send(&message, 1, MPI_INT, next,0, MPI_COMM_WORLD);
    }
    free(array), array =NULL;
    return;
}

/*A process find a solution, in total there are 12 solutions, but there is a particularity.
 * Starting from one of the solutions we can get the other 11 in this way:
 * Solution #2 -> flip solution #1
 * Solution #3-#5-#7-#9-#11  rotating previous odd solution->
 * Solution #4-#6-#8-#10-#12 flipping the previous solution
 * */
void find_all_solutions(int *array){

    n_solutions++;

    //Solution 1 : founded by the algorithm
    add_solution(array);
    //Solution 2 : is the flipped solution 1
    find_flip(array);


    //For a solution's configuration, there are 5 other solutions that can be found just rotating solution 1
    for (int rot = 1; rot <= 5; ++rot) {

        int *rotation;
        rotation = (int*) malloc((M-1) * sizeof(int));

        //the array is a representation of the hexagon, so we recompose it in order to get the "snake" that represent the order of the rotation

        //      00  01  02                                       07  03  00             (    00  01  02    )
        //    03  04  05  06                                   12  08  04  01          (   11  12  13  03   )
        //  07  08  09  10  11     ----- trasform in ---->   16  13  09  05  02       (  10  17  18  14  04  )
        //    12  13  14  15                                   17  14  10  06          (   09  16  15  05   )
        //      16  17  18                                       18  15  11             (    08  07  06    )
        transform_array_in_rotation(array,rotation);

        int *new_array;
        new_array = (int*) malloc((M-1) * sizeof(int));

        // extern hexagon 1
        // minimum is the value of the difference (end_triplet - start triplet) = 2 -0 = 2
        // maximum is the number of the elements of all the triplets = 4 * 3 = 12
        int min = 2;
        int max = 12;
        for (int i = min; i < max; ++i) {
            new_array[i] = rotation[i - 2];
        }
        for (int i = 0; i < min; ++i) {
            new_array[i] = rotation[i + max - min];
        }
        // middle hexagon 2
        // minimum is now the index of the first item remained (12 + 1 = 13)
        // maximum is the index of the last item of the middle hexagon
        min = 13;
        max = 17;

        for (int i = 12; i <= min; ++i) {
            new_array[i] = rotation[i + max - min + 1];
        }

        for (int i = min; i <= max; ++i) {
            new_array[i] = rotation[i - 1];
        }
        // inner hexagon 3
        // only 1 element
        // this element never change in every solution (invariant to rotation and flipping both)

        new_array[18] = rotation[18];

        //returning to the initial hexagonal form only for representational purposes
        transform_rotation_in_array(array, new_array);

        n_solutions ++;

        add_solution(array);
        find_flip(array);

        free(rotation), rotation = NULL;
        free(new_array), new_array = NULL;
    }

    return;
}

//this function validate the uniqueness of the values and the respect of the range constraint {1 ... 19}
bool validation(int *array){

    /*availability list*/
    bool available[19] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};

    /*any value found is not more available*/
    for (int i = 0; i < 19; ++i) {
        if(array[i] > 19 || array[i] < 1)
            return false;

        if (available[array[i]-1] == true)
            available[array[i]-1] = false;
        else return false;
    }
    return true;
}



// given in input a permutation of the 7 unknown variables, it calls back substitution on the complete matrix
void valid_permutation(int permutation[], int size){

    /* see DOCUMENTATION: improvement 1*/
    if((permutation[0] + permutation[3] + permutation [6]  > 38) ||
       (permutation[1] + permutation[3] + permutation [5]  > 38) ||
       (permutation[2] + permutation[3] + permutation [4]  > 38) )
        return;
    int k = 0;
    //restart form the initial gaussian matrix at each iteration
    int new_matrix[M-1][M];
    /*copy initial gaussian matrix in new_matrix*/
    memcpy(&new_matrix, &gaussian, ((M-1)*M)*sizeof(int));

//filling the 0-element of the main diagonal
    for (int j = 0; j < size; j++) {
        new_matrix[empty_diagonal[k]][empty_diagonal[k]] = 1;
        new_matrix[empty_diagonal[k]][M - 1] = permutation[j];
        k++;
    }
    //now call back-substitution; valid is a boolean that return constraints' respect
    bool valid = back_substitution(new_matrix);

    if (valid)
        get_final_values(new_matrix);
    return;
}

//function to swap 2 int by their address (call it swap_int(&a, &b))
void swap_int(int *a, int *b)
{
    int temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

//permutation is a recursive function
//create all the permutation of the starting_arr of int
void permutation(int *starting_arr, int start, int end)
{
    MPI_Status status;
    int received_message;

    int flag = 0; /*flag signals incoming messages if it is equal to 1*/
    /* Was chosen to call an Iprobe at every 500 iterations because of the delay caused by an Iprobe at any iterations
     * 2.38 seconds ->  2 processors, 1 Iprobe on 1 iteration
     * vs 0.196 s    -> 2 processors, 1 Iprobe on 100 iterations
     * vs 0.182 s    -> 2 processors, 1 Iprobe on 500 iterations
     * (0.300 s sequential version)*/
    if(!(iterations % 500)){
        /*check incoming messages*/
        /*previous in the ring connection, i.e. the sender*/
        int previous = (id-1) %p;
        if(id == 0)
            previous = p-1;
        MPI_Iprobe(previous,0, MPI_COMM_WORLD, &flag, &status);
        if(flag){ /*if a message was founded*/
            /*Call MPI_Recv because a message is known to be incoming*/
            MPI_Recv(&received_message,1,MPI_INT,previous,0,MPI_COMM_WORLD, &status);
            founded = 1; //because there is a message, it means that one process found the solution
            founder = received_message;
            //The receiver of the ring network connection is id+1%p
            int next = (id+1) % p;
            // but we fix the problematic case -> next = 0
            if (id == p-1)
                next = 0;
            if (id + 1 != founder) /*It is not needed to communicate again to the founder*/
                /*send in the next of ring network*/
                MPI_Send(&received_message, 1, MPI_INT, next ,0, MPI_COMM_WORLD);
        }
    }

    if (start == end) {
        if(iterations ==  analyzed + p) { // data decomposition
            analyzed = iterations;        //updating last iteration analyzed
            valid_permutation(starting_arr, end + 1);
        }
        iterations++;

        return;
    }
        for (int i = start; i <= end; i++) {
            if (!founded) {
                //swapping numbers
                swap_int((starting_arr + i), (starting_arr + start));
                /* fixing one first digit
                 * and calling permutation on
                 * the rest of the digits*/
                permutation(starting_arr, start + 1, end);
                swap_int((starting_arr + i), (starting_arr + start));
            }
        }
}

/*add a solution (array of 19 elements) in the matrix of all solutions [#solutions][19]*/
void add_solution(int *solution){

    /*if is the first solution found, allocate the memory*/
    if (n_solutions == 1){
        solutions = (int **)malloc(n_solutions * sizeof(int*));
        solutions[0] = (int *)malloc((M-1) * sizeof(int));

    }
    else{ //realloc
        solutions = realloc(solutions, (n_solutions)*sizeof(int*));
        solutions[n_solutions-1] = (int *)malloc((M-1) * sizeof(int));

    }

    for (int i = 0; i < M-1; ++i) {
        solutions[n_solutions-1][i] = solution[i];
    }




}

/*utility function that transform row-major-order representation of the hexagon in a "snake" representation"*/
void transform_array_in_rotation (int *array, int *rotation){
    rotation[0] = array[0];
    rotation[1] = array[1];
    rotation[2] = array[2];
    rotation[3] = array[6];
    rotation[4] = array[11];
    rotation[5] = array[15];
    rotation[6] = array[18];
    rotation[7] = array[17];
    rotation[8] = array[16];
    rotation[9] = array[12];
    rotation[10] = array[7];
    rotation[11] = array[3];
    rotation[12] = array[4];
    rotation[13] = array[5];
    rotation[14] = array[10];
    rotation[15] = array[14];
    rotation[16] = array[13];
    rotation[17] = array[8];
    rotation[18] = array[9];


}

/*utility function that transform "snake" representation of the hexagon in a "row-major-order" representation"*/
void transform_rotation_in_array (int *array, int *new_array) {

    array[0] = new_array[0];
    array[1] = new_array[1];
    array[2] = new_array[2];
    array[3] = new_array[11];
    array[4] = new_array[12];
    array[5] = new_array[13];
    array[6] = new_array[3];
    array[7] = new_array[10];
    array[8] = new_array[17];
    array[9] = new_array[18];
    array[10] = new_array[14];
    array[11] = new_array[4];
    array[12] = new_array[9];
    array[13] = new_array[16];
    array[14] = new_array[15];
    array[15] = new_array[5];
    array[16] = new_array[8];
    array[17] = new_array[7];
    array[18] = new_array[6];
}


/*starting from a solution, is calculated the flipped one
 *
 *           a  b  c                                     q  r  s
 *          d  e  f  g                                  m  n  o  p
 *        h  i  j  k  l     ----- trasform in ---->   h  i  j  k  l
 *          m  n  o  p                                  d  e  f  g
 *           q  r  s                                     a  b  c
 * */
void find_flip(int *array){
    n_solutions++;//obviously new solution is founded
    int *flipped_array;
    flipped_array = (int*) malloc((M-1) * sizeof(int));

    if(!flipped_array){
        fprintf(stderr,"Malloc error\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    for (int i = 0; i < M-1; ++i) {
        flipped_array[i] = array[i];
    }

    //flipping the first row to the last
    for (int i = 0; i <= 2; ++i) {
        swap_int(&flipped_array[i], &flipped_array[16+i]);

    }
    //flipping the second row to the last but one
    for (int i = 3; i <= 6; ++i) {
        swap_int(&flipped_array[i], &flipped_array[i + 9]);

    }
    add_solution(flipped_array);
    free(flipped_array), flipped_array = NULL;
}

/*utility function that read a rows x cols matrix from a file .txt*/
int read_matrix(int rows, int cols, int (*a)[cols], const char* filename)
{

    FILE *pf;
    pf = fopen (filename, "r");

    if(!pf){
        fprintf(stderr,"Error opening file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < cols; ++j)
            fscanf(pf, "%d", a[i] + j);
    }

    fclose (pf);
    return 1;
}

/*Utility function that print the final results*/
void print_hexagons(){

    for (int i = 0; i < n_solutions; ++i) {
        int *solution;
        solution= (int*) malloc((M-1) * sizeof(int));
        for (int j = 0; j < M-1; ++j) {
            solution[j] = solutions[i][j];
        }

        char solution_char[M-1][3] = {"00", "00", "00","00", "00", "00","00", "00", "00","00", "00", "00","00", "00", "00","00", "00", "00", "00"};

        for (int i = 0; i < M-1; ++i) {
            if (solution[i] < 10){ //fixing print problems
                sprintf(solution_char[i], "0%d", solution[i]);
            }
            else
                sprintf(solution_char[i], "%d", solution[i]);
        }
        fprintf(stdout,"    ***  Solution # %d  ***\n\n",i+1);

        fprintf(stdout, "          %s  %s  %s       \n", solution_char[0], solution_char[1], solution_char[2]);
        fprintf(stdout, "        %s  %s  %s  %s     \n", solution_char[3], solution_char[4], solution_char[5], solution_char[6]);
        fprintf(stdout, "      %s  %s  %s  %s  %s    \n", solution_char[7], solution_char[8], solution_char[9], solution_char[10],solution_char[11]);
        fprintf(stdout, "        %s  %s  %s  %s     \n", solution_char[12], solution_char[13], solution_char[14], solution_char[15]);
        fprintf(stdout, "          %s  %s  %s       \n\n", solution_char[16], solution_char[17], solution_char[18]);

        free(solution),solution = NULL;
    }
}

/*swap 2 rows of A matrix knowing 2 row's indexes*/
void swap_row(int idx1, int idx2)
{
    int temp;
    for (int i = 0; i < M; ++i) {
        temp = A[idx1][i];
        A[idx1][i] = A[idx2][i];
        A[idx2][i] = temp;
    }
}

//*multiply a row i of matrix A for a constant c
void multiply(int i, int c){
    for (int j = 0; j < M; ++j) {
        A[i][j] = c*A[i][j];
    }
}

/*add c*A[j][x] to each element of a row i of matrix A
 x is the index of the column*/
void add(int i, int j, int c){
    for (int x = 0; x < M; ++x) {
        A[i][x] = A[i][x] + c*A[j][x];
    }
}
