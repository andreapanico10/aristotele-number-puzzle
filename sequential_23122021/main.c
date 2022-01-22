#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

typedef enum { false, true } bool;
double time_spent = 0.0;
int n = 19;
int m = 20;
int n_solutions = 0;
bool founded = false;
int *empty_diagonal;

int gaussian[19][20];
int A[15][20];
int **solutions;

void swap(int i, int j, double a[15][20]);
void print_matrix(int matrix[n][m], int n, int m);
bool back_substitution(double matrix[n][m]);
void get_final_values (double gaussian[n][m]);
void print_results(int *array);
bool validation(int *array);
void swap_int(int *a, int *b);
void permutation(int *a, int l, int r);
void add_solution(int *solution);
void insert_permutation(int *permutation, int size);
void transform_array_in_rotation(int *array, int *rotation);
void transform_rotation_in_array(int *array, int *rotation);
void find_flip(int *array);
void getCombination(int arr[], int n,int r);
void combinationUtil(int arr[], int data[], int start, int end,int index, int r);
int readmatrix(int rows, int cols, int (*a)[cols], const char* filename);
void print_hexagons();

int main() {

    struct timeval start, end;
    int arr[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    int r = 7;
    int n = sizeof(arr)/sizeof(arr[0]);
    gettimeofday(&start, NULL);

    readmatrix(19,20,gaussian, "/Users/andrea/Desktop/gaussian_matrix.txt");
    readmatrix(15,20, A,"/Users/andrea/Desktop/initial_matrix.txt" );

    //keep track of the index of the variables unassigned in the gaussian form;
    int elements = 1;

    empty_diagonal= (int*) malloc(elements * sizeof(int));
    for (int i = 0; i < n; ++i) {
        if (gaussian[i][i] == (int) 0) {
            elements++;
            empty_diagonal = realloc(empty_diagonal, elements * sizeof(int));
            empty_diagonal[elements -2 ] = i;
        }
    }
    getCombination(arr, n, r);
    gettimeofday(&end, NULL);
    double time_taken = end.tv_sec + end.tv_usec / 1e6 - start.tv_sec - start.tv_usec / 1e6; // in seconds
    print_hexagons();
    fprintf(stdout,"Execution time: %f seconds\n", time_taken);
    free(solutions), solutions = NULL;
    return 0;
}

// function that get all combinations of size r
// in arr[] of size n. This function mainly uses combinationUtil()
void getCombination(int arr[], int n, int r)
{
    // A temporary array to store all combination one by one
    int data[r];

    // Print all combination using temporary array 'data[]'
    combinationUtil(arr, data, 0, n-1, 0, r);
}

/* arr[]  ---> Input Array
   data[] ---> Temporary array to store current combination
   start & end ---> Starting and Ending indexes in arr[]
   index  ---> Current index in data[]
   r ---> Size of a combination to be printed */
void combinationUtil(int arr[], int data[], int start, int end,int index, int r){
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        if(!founded){
            permutation(data, 0, 6);
            return;}
        else
            return;
    }

    // replace index with all possible elements. The condition
    // "end-i+1 >= r-index" makes sure that including one element
    // at index will make a combination with remaining elements
    // at remaining positions

    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1, end, index+1, r);
    }
}


void swap(int i, int j ,double A[15][20]){

    for (int k = 0; k < m; ++k) {
        double temp = A[i][k];
        A[i][k] = A[j][k];
        A[j][k] = temp;

    }


}

//utility function to print a matrix A of dimensions n*m
void print_matrix(int A[n][m], int n, int m){


    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < m; ++j)
            printf("%-3d ", A[i][j]);
        puts("");
    }
}

// this method find upper triangular matrix' n solutions
bool back_substitution(double matrix[n][m]){

    for (int i = n-1; i >= 0 ; --i) {
        double value = matrix[i][m-1]/matrix[i][i];
        //excluding all the solutions xi that have value b <= 0 || > 19
        if(value >19.0 || value < 1.0)
            return false;

        for (int j = 0; j < i; ++j) {
            if(matrix[j][i] != 0){
                matrix[j][m-1] = matrix[j][m-1] - (value * matrix[j][i]);
                matrix[j][i] = 0;
            }
        }
    }
    return true;
}

void get_final_values(double gaussian[n][m]){

    int *array= (int*) malloc(n * sizeof(int));
    int index = 0;

    // array is the list of the known terms,i.e. the last column
    for (int i = 0; i < m; ++i) {
        array[index] = (int)gaussian[i][m-1];
        index++;
    }

    //calling the validation function that check the constraints
    bool valid = validation(array);
    if(valid == true){
        //if there is no violation -> print the results
         print_results(array);
    }
    free(array), array = NULL;
    return ;
}

void print_results(int *array){
    n_solutions++;

    //Solution 1 : founded by the algorithm
    add_solution(array);
    //Solution 2 : is the flipped solution 1
    find_flip(array);

    int *rotation;
    rotation = (int*) malloc(n * sizeof(int));

    //For a solution's configuration, there are 5 other solutions that can be found just rotating solution 1
    for (int rot = 1; rot <= 5; ++rot) {
        //the array is a representation of the hexagon, so we recompose it in order to get the "snake" that represent the order of the rotation

        //      00  01  02                                       07  03  00             (    00  01  02    )
        //    03  04  05  06                                   12  08  04  01          (   11  12  13  03   )
        //  07  08  09  10  11     ----- trasform in ---->   16  13  09  05  02       (  10  17  18  14  04  )
        //    12  13  14  15                                   17  14  10  06          (   09  16  15  05   )
        //      16  17  18                                       18  15  11             (    08  07  06    )
        transform_array_in_rotation(array,rotation);

        int *new_array;
        new_array = (int*) malloc(n * sizeof(int));

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
            new_array[i] = rotation[i + max - min];
        }


        for (int i = min; i < max; ++i) {
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
        free(new_array), new_array = NULL;
    }
    founded = true;
    free(rotation), rotation = NULL;
    return;

}

//this function validate the unicity of the values and the respect of the range constraint (1 ... 19)
bool validation(int *array){

    bool available[19] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};

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
void insert_permutation(int permutation[], int size){

    int k = 0;
    double new_matrix[19][20] = {
            //a  b  c  d  e  f  g  h  i   j   k  l  m   n   o   p  q   r   s
            /*a*/{1, 0, 0, 0, 0, 0, 0, 0, 0,  1,  1, 0, 0,  1,  2,  1, 0,  1,  1,  76},
            /*b*/{0, 1, 0, 0, 0, 0, 0, 0, 0, -1,  0, 0, 0, -1, -1,  0, 0,  0,  0,  0},
            /*c*/{0, 0, 1, 0, 0, 0, 0, 0, 0,  0, -1, 0, 0,  0, -1, -1, 0, -1, -1, -38},
            /*d*/{0, 0, 0, 1, 0, 0, 0, 0, 0, -1, -1, 0, 0,  0, -1,  0, 0,  0,  0,  0},
            /*e*/{0, 0, 0, 0, 1, 0, 0, 0, 0,  0, -1, 0, 0, -1, -1, -1, 0, -1,  0, -38},
            /*f*/{0, 0, 0, 0, 0, 1, 0, 0, 0,  1,  1, 0, 0,  1,  1,  1, 0,  0,  0,  38},
            /*g*/{0, 0, 0, 0, 0, 0, 1, 0, 0,  0,  1, 0, 0,  0,  1,  0, 0,  1,  0,  38},
            /*h*/{0, 0, 0, 0, 0, 0, 0, 1, 0,  0,  0, 0, 0, -1, -1, -1, 0, -1, -1, -38},
            /*i*/{0, 0, 0, 0, 0, 0, 0, 0, 1,  1,  1, 0, 0,  1,  1,  0, 0,  1,  0,  38},
            /*j*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 0,  0,  0,  0},
            /*k*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 0,  0,  0,  0},
            /*l*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 1, 0,  0,  0,  1, 0,  0,  1,  38},
            /*m*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 1,  1,  1,  1, 0,  0,  0,  38},
            /*n*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 0,  0,  0,   0},
            /*o*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 0,  0,  0,   0},
            /*p*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 0,  0,  0,   0},
            /*q*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 1,  1,  1,   38},
            /*r*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 0,  0,  0,   0},
            /*s*/{0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, 0, 0,  0,  0,  0, 0,  0,  0,   0},
    };

//filling the 0-element of the main diagonal
    for (int j=0; j<size; j++){
        new_matrix[empty_diagonal[k]][empty_diagonal[k]] = 1;
        new_matrix[empty_diagonal[k]][m-1] = permutation[j];
        k++;
    }
    //now call back-substitution
    bool esito = back_substitution(new_matrix);

    if(esito){
        get_final_values(new_matrix);
    }
}

//function to swap 2 int by their address (call it swap_int(&a, &b))
void swap_int(int *a, int *b)
{
    int temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

//permutation recursive function
//create all the permutation of the starting_arr of int
void permutation(int *starting_arr, int start, int end)
{
    if(start==end )
    {
        if(!founded){
            insert_permutation(starting_arr, end+1);
            return;}
        else
            return;
    }
    int i;
    for(i=start;i<=end;i++)
    {
        //swapping numbers
        swap_int((starting_arr+i), (starting_arr+start));
        //fixing one first digit
        //and calling permutation on
        //the rest of the digits
        permutation(starting_arr, start+1, end);
        swap_int((starting_arr+i), (starting_arr+start));
    }
}

void add_solution(int *solution){

    if (n_solutions == 1){
        solutions = (int **)malloc(n_solutions * sizeof(int*));
        solutions[0] = (int *)malloc(n * sizeof(int));

    }
    else{
        solutions = realloc(solutions, (n_solutions)*sizeof(int*));
        solutions[n_solutions-1] = (int *)malloc(n * sizeof(int));
    }

    for (int i = 0; i < n; ++i) {
        solutions[n_solutions-1][i] = solution[i];
    }

}

void print_array(int *array, int size){
    for (int i = 0; i < size; ++i) {
        fprintf(stdout,"%d ", array[i]);

    }
    fprintf(stdout,"\n");
}

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


void find_flip(int *array){
    n_solutions++;
    int *flipped_array;
    flipped_array = (int*) malloc(n * sizeof(int));

    for (int i = 0; i < 19; ++i) {
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

int readmatrix(int rows, int cols, int (*a)[cols], const char* filename)
{

    FILE *pf;
    pf = fopen (filename, "r");
    if (pf == NULL)
        return 0;

    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < cols; ++j)
            fscanf(pf, "%d", a[i] + j);
    }

    fclose (pf);
    return 1;
}


void print_hexagons(){

    for (int i = 0; i < n_solutions; ++i) {
        int *solution;
        solution= (int*) malloc(n * sizeof(int));

        for (int j = 0; j < n; ++j) {
            solution[j] = solutions[i][j];
        }


        char solution_char[19][3] = {"00", "00", "00","00", "00", "00","00", "00", "00","00", "00", "00","00", "00", "00","00", "00", "00", "00"};

        for (int i = 0; i < 19; ++i) {
            if (solution[i] < 10){
                sprintf(solution_char[i], "0%d", solution[i]);
            }
            else
                sprintf(solution_char[i], "%d", solution[i]);
        }
        fprintf(stdout,"    ***  Solution # %d  ***\n\n",i+1);

        fprintf(stdout, "          %s  %s  %s       \n", solution_char[0], solution_char[1], solution_char[2]);
        fprintf(stdout, "        %s  %s  %s  %s     \n", solution_char[3], solution_char[4], solution_char[5], solution_char[6]);
        fprintf(stdout, "      %s %s  %s  %s  %s    \n", solution_char[7], solution_char[8], solution_char[9], solution_char[10],solution_char[11]);
        fprintf(stdout, "        %s  %s  %s  %s     \n", solution_char[12], solution_char[13], solution_char[14], solution_char[15]);
        fprintf(stdout, "          %s  %s  %s       \n\n", solution_char[16], solution_char[17], solution_char[18]);

        free(solution),solution = NULL;
     }

}