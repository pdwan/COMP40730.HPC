/* ********************************************************************************
*
*   Paula Dwan
*   13208660 : COMP40730 : Assignment 1 : Part I
*
*   STRAIGHT-FORWARD IJK IMPLEMENTATION
* 
*   Write C programs implementing the following three algorithms of multiplication of two n×n 
*   dense matrices:
*   1)    Straightforward non-blocked ijk algorithm. <*****
*   2)    Blocked ijk algorithm using square b×b blocks.
*   3)    Blocked kij algorithm using square b×b blocks.
*
*   Experiment with the programs and build/plot:
*   1)  The dependence of the execution time of each program on the matrix size n and the 
*        block size b .
*   b.  BLAS calls
*   2)  The speedup of the blocked algorithms over the non-blocked one as a function of the  
*        matrix size and the block size.
*   a.  manually written      
*   3)  Compare the fastest program with the BLAS/ATLAS routine dgemm implementing the 
*        same operation. 
*        Explain the results.
*   c.  ATLAS      
*
* *********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>

// function to obtain value from user 
int get_value (int max_value)
{
    int l_value=max_value;

    scanf("%d", &l_value);
    if (l_value >= max_value)
    {
    printf("\tWARNING : %d was entered. Now set to <max_n_value>=%d.\n", l_value, max_value);
    l_value = max_value;
    } 
    else
    {
        printf("# Entered : %d \n", l_value);
   }
   return l_value;
}

double** allocate_memory_2Dmatrix (int rows, int cols)
{
    double** local_matrix;  
    int i=0;

    local_matrix = (double**) malloc(rows * sizeof(double*));  
    if(!(local_matrix))
    {
        printf("\tERROR : exiting - memory allocation failed for matrix.\n");
        exit(1);
    }
    for (i=0; i<rows; i++)  
    {
        local_matrix[i] = (double*) malloc(cols * sizeof(double));  
        if (!(local_matrix[i]))
        {
            printf("\tERROR : exiting - memory allocation failed for matrix.\n");
            exit(1);
        }
    }
    return local_matrix;  
} 

void deallocate_2Dmatrix_memory(double **matrix2d, int rows)
{
    int ni=0;
    for (ni=0; ni<rows; ni++)
    {
        free(matrix2d[ni]);
    }
    free(matrix2d); 
}

void init_2Dmatrix_random(double **matrix2d, int rows, int cols)
{
    int ni, nj, mod_rand=10;

    printf("\n|A| or |B| : Matrix initialization using (mod %d)+1 on random number for matrix of size <%d> x <%d> ...\n", mod_rand, rows, cols);
    for(ni=0; ni<rows; ni++)
    {
        printf("[%d]\t", ni);
    }
    printf("\n");
    for(ni=0; ni<rows; ni++)
    {
        for(nj=0; nj<cols; nj++)
        {
            matrix2d[ni][nj] =  (double) ((rand()%mod_rand) + 1) ;
            printf("%g\t", matrix2d[ni][nj]); 
        }
    printf ("\n");
    }
}

void init_2Dmatrix_row(double **matrix2d, int rows, int cols)
{
    int ni, nj;
    
    printf("\n|A| or |B| : Matrix initialization with row value for matrix of size <%d> x <%d> ...\n", rows, cols);
    for(ni=0; ni<rows; ni++)
   {
        printf("[%d]\t", ni);
    }
    printf("\n");
    for(ni=0; ni<rows; ni++)
    {
        for(nj=0; nj<cols; nj++)
        {
          matrix2d[ni][nj] =  ni;
           printf("%g\t", matrix2d[ni][nj]); 
        }
        printf ("\n");
    }
}

void init_C_2Dmatrix(double **matrix2d, int rows, int cols)
{
    int ni, nj;
    
    printf("\nC : Matrix initialization with value of zero for matrix of size <%d> x <%d> ...\n", rows, cols);
    for(ni=0; ni<rows; ni++)
   {
        printf("[%d]\t", ni);
    }
    printf("\n");
    for(ni=0; ni<rows; ni++)
    {
    for(nj=0; nj<cols; nj++)
    {
        matrix2d[ni][nj] =  0.0;
        printf("%g\t", matrix2d[ni][nj]); 
    }
      printf ("\n");
   }
}

void multipy_ABC_cblas(double **matrix_a, double **matrix_b, double **matrix_c, int rows, int cols, double alpha, double beta)
{
// m, n, k :   local integers indicating the size of the matrices:
//                 A :  m rows by k columns
//                 B :  k rows by n columns
//                 C :  m rows by n columns
//                 Here, m = n = k = rows = columns = <nx> as entered
   int lm = rows, ln = rows;
   int ni, nj;
// la_offset, lb_offset, lc_offset :
//                  Leading dimension of matrix A, B or C respectively, or the number of elements between 
//                  successive rows for row-major storage or columns for column-major storage. 
   int la_offset = rows, lb_offset = cols, lc_offset = rows;

   printf("   Simple ijk for C[i][j] using cblas ... \n");
   # cblas_dgemm( CblasRowMajor,  CblasTrans, CblasNoTrans, lm, ln, ln, alpha, matrix_a, la_offset, matrix_c, lb_offset, beta, matrix_c, lc_offset);
   
   printf("   Calculation completed - results as follows ... \n");

// NB compare like with like calculation and output : write to file or stdout for both manual and dgemm
   for(ni=0; ni<rows; ni++)
   {
      printf("[%d]\t", ni);
   }
   for (ni=0; ni<rows; ni++)
   {
      printf("   ");
      for (nj=0; nj<cols; nj++)
      {
         printf("%g\t",matrix_c[ni][nj]); 
      }
   printf("\n");
   }
}

void multipy_ABC_manual(double **matrix_a, double **matrix_b, double **matrix_c, int rows, int cols)
{
   int ni, nj, nk;

   for(ni=0; ni<rows; ni++)
   {
      printf("[%d]\t", ni);
   }
   printf("\n");
   for (ni=0; ni<rows; ni++)
   {
      printf("   ");
      for (nj=0; nj<cols; nj++)
      {
         double sum = 0.0;
         for (nk=0; nk<rows; nk++)
         {
            sum+= (matrix_a[ni][nk]) * (matrix_b[nk][nj]);
         }
         matrix_c[ni][nj] = sum;
         printf("%g\t",matrix_c[ni][nj]); 
      }
   printf("\n");
   }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main (void)
{
// define variables
   struct timeval tv1, tv2;
   struct timezone tz;
   const int MAXNVALUE = 100;
   const double ALPHA = 1.0;
   const double BETA = 1.0;

// start of program
   printf("\n......................................................................................................... \n");
   printf("INPUT :\tstraight-forward ijk (here nx = ny, i.e. |row| = |col|) ...\n");  
   printf("\tPlease enter matrix dimension less than <max_value> = %d --> <nx> : ", MAXNVALUE);
   int nx = get_value(MAXNVALUE);
   int ny = nx; 

// Allocate memory for matrices A, B & C
   printf("CREATE MATRICES : \n");
   double **A =  allocate_memory_2Dmatrix(nx, ny);
   double **B =  allocate_memory_2Dmatrix(nx, ny);
   double **C =  allocate_memory_2Dmatrix(nx, ny);

// Initialize matrices A & B
   printf("INITIALIZE MATRICES : \n");
   init_2Dmatrix_random(A, nx, ny);
   init_2Dmatrix_random(B, nx, ny);
   //init_2Dmatrix_row(A, nx, ny);
   //init_2Dmatrix_row(B, nx, ny);
   (C, nx, ny);

// Calculate program execution time taken : start
   gettimeofday(&tv1, &tz);
// Execute multiplication of straightforward ijk 2-D matrix NxN
   printf("RESULTS : manual calculation - \n");
   multipy_ABC_manual(A, B, C, nx, ny);
// program execution time taken : end
   gettimeofday(&tv2, &tz);
   double elapsed_manual = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

//  Calculate program execution time taken : start
    gettimeofday(&tv1, &tz);
//  Execute multiplication of straightforward ijk 2-D matrix NxN
    printf("RESULTS : BLAS/ATLAS calculation -\n");
    multipy_ABC_cblas(A, B, C, nx, ny, ALPHA,BETA);
//  program execution time taken : end
    gettimeofday(&tv2, &tz);
    double elapsed_cblas = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

//  output results to .dat file : matrix size || Time/manual 	|| Time / dgemm
    printf("   %g \t %d \n",  nx, elapsed_manual, elapsed_cblas);
   	
//  De-allocate memory for matrices A, B & C
    printf("CLEAN-UP : \n");
    deallocate_2Dmatrix_memory(A, nx);
    deallocate_2Dmatrix_memory(B, nx);
    deallocate_2Dmatrix_memory(C, nx);
  
//  end of program
   printf("\n......................................................................................................... \n");

   return 0;
}

