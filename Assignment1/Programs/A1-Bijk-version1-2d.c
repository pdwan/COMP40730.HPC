/*
*
* Paula Dwan
* 13208660 : COMP40730 : Assignment 1 : Part II
*
* BLOCKED IJK IMPLEMENTATION
* 
* Write C programs implementing the following three algorithms of multiplication of two n×n dense 
* matrices:
*       1)    Straightforward non-blocked ijk algorithm.
* --> 2)    Blocked ijk algorithm using square b×b blocks.
*       3)    Blocked kij algorithm using square b×b blocks.
*
* Experiment with the programs and build/plot:
*       1)
*       The dependence of the execution time of each program on the matrix size n and the block size b .
*       b.     BLAS calls
*       2)
*       The speedup of the blocked algorithms over the non-blocked one as a function of the matrix size and 
*       the block size.
*       a.     manually written      
*       3)
*       Compare the fastest program with the BLAS/ATLAS routine  dgemm implementing the same 
*       operation. 
*       Explain the results.
*       b.     ATLAS      
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

//function to obtain value from user 
int get_value (int max_value)
{
   int l_value=max_value;

   scanf("%d", &l_value);
   if (l_value >= max_value)
   {
      printf("   Warning : <value> = %d was entered. Now set to <max_n_value>=%d.\n", l_value, max_value);
      l_value = max_value;
   } 
   else
   {
      printf("   <value> = %d was entered.\n", l_value);
   }
   return l_value;
}

// function to allocate matrix memory 
double** create_2Dmatrix (int rows, int cols)
{
   double** local_matrix;  
   int i=0;

   local_matrix = (double**) malloc(rows * sizeof(double*));  
   if(!(local_matrix))
   {
      printf("   Error : exiting - memory allocation failed for matrix.\n");
      exit(1);
   }
   for (i=0; i<rows; i++)  
   {
      local_matrix[i] = (double*) malloc(cols * sizeof(double));  
      if (!(local_matrix[i]))
      {
         printf("   Error : exiting - memory allocation failed for matrix.\n");
         exit(1);
      }
   }
   printf("   Memory allocation complete for matrix ...\n");
   return local_matrix;  
} 

// function to deallocate matrix memory
void deallocate_2Dmatrix_memory(double **arr, int rows)
{
   int ni=0;
 
   for (ni=0; ni<rows; ni++)
   {
      free(arr[ni]);
   }
   free(arr); 
   printf("   Memory de-allocation completed for matrix ...\n");
}

// function to initialize the matrices A and B using random numbers
void init_2Dmatrix_random(double **arr, int rows, int cols)
{
   int ni, nj, mod_rand=10;

   printf("   Matrix initialization using (mod %d)+1 on random number for matrix of size <%d> x <%d> ...\n", mod_rand, rows, cols);
   printf("   ");
   for(ni=0; ni<rows; ni++)
   {
      printf("[%d]\t", ni);
   }
   printf("\n");
   for(ni=0; ni<rows; ni++)
   {
      printf("   ");
      for(nj=0; nj<cols; nj++)
      {
         arr[ni][nj] =  (double) ((rand()%mod_rand) + 1) ;
         printf("%g\t", arr[ni][nj]); 
      }
      printf ("\n");
   }
}

// function to initialize the matrices A and B using specific value
void init_2Dmatrix_specific(double **arr, int rows, int cols, double element_value)
{
   int ni, nj;

   printf("   Matrix initialization with value of %g for matrix of size <%d> x <%d> ...\n", element_value, rows, cols);
   printf("   ");
   for(ni=0; ni<rows; ni++)
   {
      printf("[%d]\t", ni);
   }
   printf("\n");
   for(ni=0; ni<rows; ni++)
   {
      printf("   ");
      for(nj=0; nj<cols; nj++)
      {
         arr[ni][nj] =  element_value;
         printf("%g\t", arr[ni][nj]); 
      }
      if (element_value != 0.0) 
      {
         element_value ++;
      }
      printf ("\n");
   }
}

//function to multiply 2-D matrices A, B and C using Blocked ijk
void multipy_ABC_coded_simple(double **arr_a, double **arr_b, double **arr_c, int rows, int cols, int block)
{
   int ni, nj, nk;

   printf("   Blocked (simplistic calculation : <nx> = <ny> div <nb>) ijk for C[ni][nj] ... \n");

   for (ni=0; ni<(rows/block); ni++)
   {
      printf("   ");  
      for (nj=0; nj<(rows/block); nj++) 
      {
         printf("C[%d][%d]\t", ni, nj);
      }
      printf("\n");
   }
   for (ni=0; ni<(rows/block); ni++)
   {
      printf("   ");  
      for (nj=0; nj<(rows/block); nj++) 
      {
         double sum = 0.0;
         for (nk=0; nk<(cols/block); nk++)
         {
            sum += arr_a[ni][nk] * arr_b[nk][nj];
         }
         arr_c[ni][nj] = sum;
         printf("%g\t", arr_c[ni][nj]);
      }
      printf("\n");
   }
}


//function to multiply 2-D matrices A, B and C using Blocked ijk
void multipy_ABC_coded(double **arr_a, double **arr_b, double **arr_c, int rows, int cols, int block)
{
   int bi, bj, bk, ni, nj, nk;

   printf("   Blocked ijk for C[bi+ni][bj+nj] ... \n");
   for (bi=0; bi<rows; bi+=block)
   {
      for (bj=0; bj<cols; bj+=block) 
      {
         for (bk=0; bk<rows; bk+=block) 
         {
            for (ni=0; ni<block; ni++)
            {
               printf("   ");  
               for (nj=0; nj<block; nj++) 
               {
                  double sum = 0.0;
                  for (nk=0; nk<block; nk++)
                  {
                     sum += arr_a[bi+ni][bk+nk] * arr_b[bk+nk][bj+nj];
                  }
                  arr_c[bi+ni][bj+nj] = sum;
                  printf("C[%d+%d][%d+%d] = %g\t", bi, ni, bj, nj, sum);
               }
               printf("\n");
            }
         }
      }
   }
}

//function to multiply 2-D matrices A, B and C using Blocked ijk
void multipy_ABC_cblas(double **arr_a, double **arr_b, double **arr_c, int rows, int cols, int block)
{

  printf("   Blocked ijk for C[bi+ni][bj+nj] using cblas ... \n");

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//start of main

int main (void)
{
// define variables
   struct timeval tv1, tv2;
   struct timezone tz;
   const int MAXNVALUE = 100;

// start of program
   printf("\n");
   printf("\n................................................................................................ \n");
   printf("\nSTART :\n");

   printf("INPUT : \n   Blocked ijk (here nx = ny, i.e. |row| = |col|) ...\n");

// Get n value for row x col matrix, both equal to n
   printf("   Please enter matrix dimension less than <max_value> = %d --> <nx> : ", MAXNVALUE);
   int nx = get_value(MAXNVALUE);
   int ny = nx;

// Get block size
   printf("   Please enter sub-matrix dimension i.e.: block size, less than <max_value> = %d --> <nb> : ", MAXNVALUE);
   int nb = get_value(MAXNVALUE);
   while ( (!((nx % nb)==0)) && (!(nb<nx)) )
   {
      printf("   Error : (i) nb not < nx and (nx mod nb) is not equal to 0. Please re-enter <nb> : ");
      nb = get_value(MAXNVALUE);
   }

// Allocate memory for matrices A, B & C
   printf("\nCREATE MATRICES : \n");
   double **A = create_2Dmatrix(nx, ny);
   double **B = create_2Dmatrix(nx, ny);
   double **C = create_2Dmatrix(nx, ny);

// Initialize matrices A & B
   printf("\nINITIALIZE MATRICES : \n");
   init_2Dmatrix_random(A, nx, ny);
   init_2Dmatrix_random(B, nx, ny);
   //init_2Dmatrix_specific(A, nx, ny, 2.0);
   //init_2Dmatrix_specific(B, nx, ny, 3.0);
   init_2Dmatrix_specific(C, nx, ny, 0.0);

// Calculate program execution time taken : start
   gettimeofday(&tv1, &tz);

// Execute multiplication of blocked ijk 2-D matrix NxN
   printf("\nRESULTS : for loops \n");
   multipy_ABC_coded(A, B, C, nx, ny, nb);

// program execution time taken : end
   gettimeofday(&tv2, &tz);
   double elapsed_manual = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

// output results, continued for time 
   printf("   Time taken to complete : %g seconds.\n", elapsed_manual);

// Calculate program execution time taken : start
   gettimeofday(&tv1, &tz);
   
// Calculate program execution time taken : start
   gettimeofday(&tv1, &tz);

// Execute multiplication of blocked ijk 2-D matrix NxN
   printf("\nRESULTS : for loop simple  \n");
   multipy_ABC_coded_simple(A, B, C, nx, ny, nb);

// program execution time taken : end
   gettimeofday(&tv2, &tz);
   double elapsed_manual_simple = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

// output results, continued for time 
   printf("   Time taken to complete : %g seconds.\n", elapsed_manual_simple);

// Calculate program execution time taken : start
   gettimeofday(&tv1, &tz);

// Execute multiplication of blocked ijk 2-D matrix NxN
   printf("RESULTS : cblas\n");
   multipy_ABC_cblas(A, B, C, nx, ny, nb);
  
// program execution time taken : end
   gettimeofday(&tv2, &tz);
   double elapsed_cblas = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

// output results, continued for time 
   printf("   Time taken to complete : %g seconds.\n", elapsed_cblas);
	
// De-allocate memory for matrices A, B & C
   printf("CLEAN-UP : \n");
   deallocate_2Dmatrix_memory(A, nx);
   deallocate_2Dmatrix_memory(B, nx);
   deallocate_2Dmatrix_memory(C, nx);
  
// end of program
   printf("\nEND : \n");
   printf("\n................................................................................................ \n");
   printf("\n");

// finish
   return 0;
}

