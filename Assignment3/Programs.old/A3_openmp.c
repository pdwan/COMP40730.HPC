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

// function to obtain value from user for matrix size, number of processes.
int get_value (int max_value)
{
     int l_value=max_value;

     scanf("%d", &l_value);
     if (l_value >= max_value)
     {
           printf("      Warning : value=%d was entered > <max_n_value>=%d, now reset to <max_n_value>. Calculation will now commence.\n", l_value, max_value);
           l_value = max_value;
     } 
     else 
     {
           printf("      n=%d was entered. Calculation will now commence using absolute value.\n", l_value);
     }
     return abs(l_value);}

// function to allocate matrix memory for 1D matrix using size of equivalent 2D matrix
double* create_1Dmatrix (int rc)
{
     double *local_matrix; 

     local_matrix = malloc(rc * sizeof(double)); 
     if(!(local_matrix))
     {
           printf("      Error : exiting - memory allocation failed for 1D matrix.\n");
           exit(1);
     }
     printf("      Memory allocation complete for 1D matrix matrix ...\n");
     return local_matrix;     
} 

// function to deallocate matrix memory
void deallocate_matrix_memory(double *arr)
{
     free(arr); 
     printf("      Memory de-allocation completed for matrix ...\n");
}

// function to initialize the non-zero 1D Matrices using random number
void init_1Dmatrix(double *arr, int rc)
{
     int ni, mod_rand=10;

     printf("      1D matrix initialization using (mod %d)+1 on random ...\n", mod_rand);
     for(ni=0; ni<rc; ni++)
     {
           arr[ni]  = (double)((rand()%mod_rand) + 1);
           printf("m[%d] = %g\t", ni, arr[ni]); 
     }
           printf ("\n");
}

// function to initialize the non-zero 1D Matrices using set value
void init_1Dmatrix_specific_value(double *arr, int rc, double element_value)
{
     int ni;
     
     printf("      1D matrix initialization using set value : %g ...\n", element_value);
     for(ni=0; ni<rc; ni++)
     {
           arr[ni] = element_value;
           printf("m[%d] = %g\t", ni, arr[ni]); 
     }
           printf ("\n");
}

// function to calculate dot product of matrix A and matrix B using slices of A
void *norm_summation_calculation(void *args)
{
     norm_summation_t *norm_summation_data; 
     int i; 
      
     printf("      Calculate dot product of <allB> and <allA> using <sliceA> writing to <sliceC>...\n");
     norm_summation_data = args; 
     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, norm_summation_data->my_nx, norm_summation_data->my_nx, norm_summation_data->my_nx, 1.0, norm_summation_data->my_sliceA, norm_summation_data->my_nx, norm_summation_data->my_allB, norm_summation_data->my_nx, 0.0, norm_summation_data->my_sliceC, norm_summation_data->my_nx);
     
     // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dot_data->my_nx, dot_data->my_nx, dot_data->my_nx, 1.0, dot_data->my_sliceA, dot_data->my_nx, dot_data->my_allB, dot_data->my_nx, 0.0, &dot_data->my_dot_prod[(i * dot_data->my_nx)], dot_data->my_nx);
          
     printf("      Evaluate sum of local <sliceC> using <mutex>...\n");
     pthread_mutex_lock(norm_summation_data->mutex);
     printf("      ................... > norm_summation_data->my_ns_size = %d \n ", (int)norm_summation_data->my_ns_size);
     printf("value \t\t  max");
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
          printf("%g \t\t %g \n",  (double)norm_summation_data->my_norm_summation[i], norm_summation_data->my_norm_sumC);
     }  
     printf("      Max absolute value in local <sliceC> is : %g\n", norm_summation_data->my_norm_sumC);     
     pthread_mutex_unlock(norm_summation_data->mutex); 
     pthread_exit(NULL); 
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// start of main
int main()
{
// define variables
     struct timeval tv1, tv2;
     struct timezone tz;
     const int MAXNMATRIX  = 100;
     const int MAXPTHREADS = 124;
     int i = 0;
          
// start of program
     printf("\n\n");
     printf("\n............................................................................................................................\n");
     printf("\nSTART :\n");

 printf("INPUT : \n      Use <pthreads> to calculate dot product and Infinity Norm of product of two NxN matrices ...\n");

// Get values for matrix row x col matrix, both equal to n
     printf("      Please enter matrix dimension < maximun of [%d] i.e.: n : ", MAXNMATRIX);
     int nx = get_value(MAXNMATRIX);
     int ny = nx;

// Get number of processes
     printf("      Please enter number of processors (i.e.: number of threads) < maximum of [%d] --> <np>: ", MAXPTHREADS);
     int np = get_value(MAXPTHREADS);
          while ( (nx % np) > 0 )
     {
           printf("      Error : matrix size = %d must  be an equal multiple of number of processors and < [%d]. Please reenter [np]: ", nx, MAXPTHREADS);
           np = get_value(MAXPTHREADS);
     }
 
// Calculate slice information using number of processors <np> and matrix size <nx>
     int ns_size = (nx*np);
     int ns_no = (nx / np); 
     printf("      No of slices = <matrix rows> x <number of processes> = [%d] x [%d] = [%d] ...\n", nx, np, ns_size);
     printf("      Slice size       = <matrix rows> ÷ <number of processes> = [%d] ÷ [%d] = [%d] ...\n", nx, np, ns_no);
 
// Allocate memory for & initialize matrices allA, allB ; create & initialize array C (row totals) 
     printf("CREATE MATRICES <allA>, <allB>, <allC>, <sliceC> & <sliceC> :\n");
     double *allA = create_1Dmatrix(nx*ny);
     double *allB = create_1Dmatrix(nx*ny);
     double *allC = create_1Dmatrix(nx*ny);
     double *sliceA = create_1Dmatrix(ns_size);
     double *sliceC = create_1Dmatrix(ns_size);
 
     printf("INITIALIZE MATRICES <allA>, <allB>, <allC>, <sliceB>  & <sliceC> :\n");
     init_1Dmatrix(allA, nx*ny);
     init_1Dmatrix(allB, nx*ny);
     init_1Dmatrix_specific_value(allC, nx*ny, 0.0);
     init_1Dmatrix_specific_value(sliceA, ns_size, 0.0);
     init_1Dmatrix_specific_value(sliceC, ns_size, 0.0);
     
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
          printf("      Create working thread[%d] and send to norm_summation_calculation function ...\n", i);
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
          printf("      Implement pthread_join for working thread[%d] ... \n", i);
          pthread_join(working_thread[i], &status); 
     }
     printf("\n      Calculate infintity norm of matrix C containing dot product of A and B using <mutex>...\n");
      
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

     printf("\nEND :\n");
     printf("\n............................................................................................................................\n");
     printf("\n\n");
  
// finish
  return 0;
}

