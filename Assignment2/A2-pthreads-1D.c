/* 
* Paula Dwan
* 13208660 : COMP40730 : Assignment 2
*
* Parallel pthreads program computing norm of product of two nxn matrices
* 1/     Granularity of the product : 
*           (ii)      One step algorithim. No intermediate resulting matrix
* 2/     Partitioning Scheme :
*           (a)      Left matrix is horizonitally partitioned
* 3/     Matrix norm to be computed :
*           (b)      Maximum absolute row sum norm (aka infinity-norm)
*                     calculate max value of sum of each row (or product of same) 
*                     ||A|| infintity = max SUM | Aij |
*                     where   max     :==> 0 <= i <=n
*                                  SUM   :==>  j =0 to n-1      
*
*/          

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <cblas.h>
#include <math.h>

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
typedef struct 
{ 
     double *my_sliceA; 
     double *my_sliceC; 
     double *my_allB; 
     double *my_norm_summation; 
     double my_norm_sumC; 
     pthread_mutex_t *mutex; 
     int my_ns_size;
     int my_nx; 
     int my_np;
} norm_summation_t;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void usage () {
    printf( "USAGE : \t<program name> <i> |<r> <N> <matrix contents file>.txt <timing file>.dat \n");
    printf ("TO : \tCalculate infinity norm of and populate matrix |C| = |A| x |B| using straight-forward IJK algorithm and cblas  - via bash script. \n");    
    printf( "where \n");
    printf("\t1. \t<R>\tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n");
    printf("\t2. \t<I>\tinitialize |A| & |B| _incrementally_ with <row> value and |C| with '0' \n");
    printf("\t3. \t<N> \tmax size of each matrix, defaults to 1,000 if invalid or not provided - verified in calling script \n");
    printf("\t4. \t<matrix contents file>.txt\n\t\tname of .txt file to store values of matrices |A| |B| & |C| \n");
    printf("\t5. \t<timing .dat file> .dat \n\t\tname of .dat file to contain time to complete for each iteration \n");
    exit(1);
}

// function to allocate matrix memory for 1D matrix using size of equivalent 2D matrix
double* create_1Dmatrix (int rc)
{
     double *local_matrix; 

     local_matrix = malloc(rc * sizeof(double)); 
     if(!(local_matrix))
     {
           printf("ERROR : \tExiting - memory allocation failed for 1D matrix.\n");
           exit(1);
     }
     printf("\tMemory allocation complete for 1D matrix matrix ...\n");
     return local_matrix;     
} 

// function to deallocate matrix memory
void deallocate_matrix_memory(double *arr)
{
     free(arr); 
     printf("\tMemory de-allocation completed for matrix ...\n");
}

void init_1d_matrix_random(double *matrix1d, int rows, int cols, FILE *fp)
{
    int ni, mod_rand=10;
    int count =1;
    fprintf(fp, "\n# |A| or |B| : \tMatrix initialization using (mod %d)+1 on random number for matrix of size <%d> x <%d> ...\n", mod_rand, rows, cols);
    for(ni=0; ni<(rows*cols); ni++)
    {
        matrix1d[ni]= (double) ((rand()%mod_rand) + 1);
        fprintf(fp, "%g\t", matrix1d[ni]); 
        if (count == rows) {
            fprintf(fp, "\n");
            count =1;
        } else 
        {
            count ++;    
        }
    }
}

void init_1d_matrix_increment(double *matrix1d, int rows, int cols, FILE *fp)
{
    int ni;    
    int count =1;
    fprintf(fp, "\n# |A| or |B| : \tMatrix initialization with row value for matrix of size <%d> x <%d> ...\n", rows, cols);
    for(ni=0; ni<(rows*cols); ni++)
    {
        matrix1d[ni]= ((rows*cols) % rows);
        fprintf(fp, "%g\t", matrix1d[ni]); 
        if (count == rows) {
            fprintf(fp, "\n");
            count =1;
        }else 
        {
            count ++;    
        }
    }
}

void init_C_1d_matrix(double *matrix1d, int rows, int cols, FILE *fp)
{
    int ni;
    int count =1;
    fprintf(fp, "\n# |C| : \tMatrix initialization with value of zero for matrix of size <%d> x <%d> ...\n", rows, cols);
    for(ni=0; ni<(rows*cols); ni++)
    {
        matrix1d[ni]=  0.0;
        fprintf(fp, "%g\t", matrix1d[ni]); 
        if (count == rows) {
            fprintf(fp, "\n");
            count =1;
        }else {
        count++;
        }
    }
}

//  function to calculate infinity norm and provide values for |C|
double multipy_ABC_complex(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols, FILE *fp)
{
    int ni, nj, nk;
    int count = 1;
    double max_row = 0.0;
    fprintf(fp, "\n# |C| : \tMatrix computed values for matrix of size <%d> x <%d> ... using SIMPLE : Blocked IJK\n", rows, cols);
    for (ni=0; ni<rows; ni++)
    {
        for (nj=0; nj<cols; nj++)
        {
            double sum = 0.0;
            for (nk=0; nk<rows; nk++)
            {
                sum+= (matrix_a[(ni*rows)+nk]) * (matrix_b[(nk*rows)+nj]);
            }
            matrix_c[ni+(nj*cols)] = sum;
            max_row = (max_row < sum) ? sum : max_row ; 
            fprintf(fp, "%g\t",matrix_c[ni+(nj*cols)]); 
            if (count == rows) 
            {
                fprintf(fp, "\n");
                count =1;
            }else 
            {
                count++;
            }
        }
    }
    return max_row;
}

//  function to calculate infinity norm and provide values for |C|
void *norm_summation_calculation(void *args FILE *fp)
{
     norm_summation_t *norm_summation_data; 
     int i;  
     printf(fp, "# Calculate infinity norm of <allB> and <allA> using <sliceA> writing to <sliceC>...\n");
     norm_summation_data = args; 
     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, norm_summation_data->my_nx, norm_summation_data->my_nx, norm_summation_data->my_nx, 1.0, norm_summation_data->my_sliceA, norm_summation_data->my_nx, norm_summation_data->my_allB, norm_summation_data->my_nx, 0.0, norm_summation_data->my_sliceC, norm_summation_data->my_nx);
     printf(fp, "# Evaluate sum of local <sliceC> using <mutex>...\n");
     pthread_mutex_lock(norm_summation_data->mutex);
     printf(fp, "# ................... > norm_summation_data->my_ns_size = %d \n ", (int)norm_summation_data->my_ns_size);
     printf(fp, "# value \t\t  max");
     int iteration = (int) (norm_summation_data->my_ns_size / norm_summation_data->my_nx);
     int counter = 0;
     for(i=0; i<iteration; i++) 
     {
          norm_summation_data->my_norm_summation[counter] += norm_summation_data->my_sliceC[i];
          if ((counter % (int) norm_summation_data->my_nx) > 0)
          { 
               counter ++;
          } else 
          {
               counter =0;
          }
          if (norm_summation_data->my_norm_sumC < norm_summation_data->my_norm_summation[i]) 
          {
               norm_summation_data->my_norm_sumC = norm_summation_data->my_norm_summation[i] ;
          }
          printf(fp, "# %g \t\t %g \n",  (double)norm_summation_data->my_norm_summation[i], norm_summation_data->my_norm_sumC);
     }  
     printf(fp, "# Max absolute value in local <sliceC> is : %g\n", norm_summation_data->my_norm_sumC);     
     pthread_mutex_unlock(norm_summation_data->mutex); 
     pthread_exit(NULL); 
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// start of main
int main ( int argc, char *argv[] )
{
//  define variables
    struct timeval tv1, tv2;
    struct timezone tz;
    const int MAXPTHREADS = 124;
    int i = 0;
    const double ALPHA = 1.0;
    const double BETA = 1.0;
    int increment_or_random = 5;
    int max_num_args = 6;
    char filename_matrix[50];
    char filename_timing[50];
              
//  start of program

    if ( argc != max_num_args )
    {
        usage();
    } 
    
    char cell_value_type[3];
    strncpy(cell_value_type, argv[1], 2);
    cell_value_type[3] = '\0';
    if (strcmp(cell_value_type,"-i") == 0)
    {
        increment_or_random = 1;        
    } else if (strcmp(cell_value_type,"-r") == 0)
    {
        increment_or_random = 0;
    }
    int nx = atoi(argv[2]);                                     
    int ny = nx;      
    int num_threads = atoi(argv[3]);                                     
    strncpy(filename_matrix, argv[4], 49);
    filename_matrix[50] = '\0';
    FILE *fp_matrix = fopen( filename_matrix, "wa" );
    strncpy(filename_timing, argv[5], 49);
    filename_timing[50] = '\0';
    FILE *fp_timing = fopen(filename_timing, "wa" );    
    fprintf(fp_matrix, "RUNNING : \t%s %s %d %d %s %s \n", argv[0],cell_value_type,nx,nb,filename_matrix,filename_timing );
    fprintf(fp_timing, "RUNNING : \t%s %s %d %d %s %s \n", argv[0],cell_value_type,nx,nb,filename_matrix,filename_timing );
 
    printf("\n............................................................................................................................\n");
    printf("INPUT :\tUse <pthreads> to calculate dot product and Infinity Norm of product of two NxN matrices ...\n");
    printf("\tPlease enter matrix dimension < maximun of [%d] i.e.: n : ", MAXNMATRIX);
    int nx = get_value(MAXNMATRIX);
    int ny = nx;
    printf("\tPlease enter number of processors (i.e.: number of threads) < maximum of [%d] --> <np>: ", MAXPTHREADS);
    int np = get_value(MAXPTHREADS);
    while ( (nx % np) > 0 )
    {
        printf("ERROR : \tmatrix size = %d must  be an equal multiple of number of processors and < [%d]. Please reenter [np]: ", nx, MAXPTHREADS);
        np = get_value(MAXPTHREADS);
    }
 
//  Calculate slice information using number of processors <np> and matrix size <nx>
     int ns_size = (nx*np);
     int ns_no = (nx / np); 
     printf("\tNo of slices = <matrix rows> x <number of processes> = [%d] x [%d] = [%d] ...\n", nx, np, ns_size);
     printf("\tSlice size       = <matrix rows> ÷ <number of processes> = [%d] ÷ [%d] = [%d] ...\n", nx, np, ns_no);
 
// Allocate memory for & initialize matrices allA, allB ; create & initialize array C (row totals) 
     printf("CREATE MATRICES <allA>, <allB>, <allC>, <sliceA> & <sliceC> :\n");
     double *allA = create_1Dmatrix(nx*ny);
     double *allB = create_1Dmatrix(nx*ny);
     double *allC = create_1Dmatrix(nx*ny);
     double *sliceA = create_1Dmatrix(ns_size);
     double *sliceC = create_1Dmatrix(ns_size);
 
 //  Initialize matrices A & B & C
    fprintf(fp_matrix, "INITIALIZE MATRICES <allA>, <allB>, <allC>, <sliceA>  & <sliceC> :\n");
    if (increment_or_random == 0) {
        init_1d_matrix_random(allA, nx, ny, fp_matrix);
        init_1d_matrix_random(allB, nx, ny, fp_matrix);
    } else if (increment_or_random == 1) 
    {
        init_1d_matrix_increment(allA, nx, ny, fp_matrix);
        init_1d_matrix_increment(allB, nx, ny, fp_matrix);
    }
    init_C_1d_matrix(allC, nx, ny, fp_matrix);
    init_C_1d_matrix(sliceA, nx, ny, fp_matrix);
    init_C_1d_matrix(sliceC, nx, ny, fp_matrix);
     
// Calculate program execution time taken : start
     gettimeofday(&tv1, &tz);

// Commence pthread calculation
     printf("RESULTS :\n");
     
     pthread_t *working_thread malloc(np * (sizeof(pthread_t)) );  
     pthread_mutex_t *mutex_dot_product = malloc(sizeof(pthread_mutex_t)); 
     norm_summation_t *thread_norm_summation_data = malloc(np * (sizeof(norm_summation_t))); 
     double norm_summation_matrix = malloc( * (sizeof(double)));
     void *status; 
     
    pthread_mutex_init(mutex_dot_product, NULL); 

    for(i=0; i<np; i++) 
    { 
        printf("\tCreate working thread[%d] and send to norm_summation_calculation function ...\n", i);
        thread_norm_summation_data[i].my_sliceA = allA + (i * ns_size);
        thread_norm_summation_data[i].my_sliceC =  allC + (i * ns_size);
        thread_norm_summation_data[i].my_allB = allB;
        thread_norm_summation_data[i].my_norm_summation = norm_summation_matrix;
        thread_norm_summation_data[i].my_norm_sumC = -1.0;
        thread_norm_summation_data[i].mutex = mutex_dot_product;
        thread_norm_summation_data[i].my_ns_size = ns_size;
        thread_norm_summation_data[i].my_nx = nx;
        thread_norm_summation_data[i].my_np = np;       
        double *my_sliceA; 
        double *my_sliceC; 
        double *my_allB; 
        double *my_norm_summation; 
        double my_norm_sumC; 
        pthread_mutex_t *mutex; 
        int my_ns_size;
        int my_nx; 
        int my_np;
        pthread_create(&working_thread[i], NULL, norm_summation_calculation,(void*)&thread_norm_summation_data[i]);
    } 
     
    for(i=0; i<np; i++) 
    {
        printf("\tImplement pthread_join for working thread[%d] ... \n", i);
        pthread_join(working_thread[i], &status); 
    }
    printf("\tCalculate infintity norm of matrix C containing dot product of A and B using <mutex>...\n");
      
// program execution time taken : end
     gettimeofday(&tv2, &tz);
     double elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

// output results, continued
     printf("      Time taken to complete calculations : %g seconds.\n", elapsed);

// De-allocate memory for matrices allA, allB & C
     printf("CLEAN-UP :\n");
     deallocate_matrix_memory(allA);
     deallocate_matrix_memory(allB);
     deallocate_matrix_memory(sliceA);
     deallocate_matrix_memory(sliceC);
     free(working_thread);
     free(thread_norm_summation_data);
     pthread_mutex_destroy(mutex_dot_product);
     free(mutex_dot_product);

     printf("\n............................................................................................................................\n");
  
// finish
  return 0;
}

