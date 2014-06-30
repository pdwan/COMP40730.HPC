/*
*
* Paula Dwan
* 13208660 : COMP40730 : Assignment 1 : Part III
*
* BLOCKED KIJ IMPLEMENTATION
* 
* Write C programs implementing the following three algorithms of multiplication of two n×n dense 
* matrices:
*       1)    Straightforward non-blocked ijk algorithm.
*       2)    Blocked ijk algorithm using square b×b blocks.
* --> 3)    Blocked kij algorithm using square b×b blocks.
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
     printf("   Warning : <value>=%d was entered. Now set to <max_n_value>=%d. Calculation will now commence.\n", l_value, max_value);
     l_value = max_value;
  } 
  else
  {
     printf("   <value>=%d was entered. Calculation will now commence.\n", l_value);
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
  for(ni=0; ni<rows; ni++)
  {
     printf("   ");
     for(nj=0; nj<cols; nj++)
     {
        arr[ni][nj] =  (double)((rand()%mod_rand) + 1);
        printf("a[%d][%d] = %g\t", ni, nj, arr[ni][nj]; 
     }
     printf ("\n");
  }
}

// function to initialize the matrices A and B using specific value
void init_2Dmatrix_specific(double **arr, int rows, int cols, double element_value)
{
  int ni, nj, mod_rand=10;

  printf("   Matrix initialization with value of %d for matrix of size <%d> x <%d> ...\n", element_value, rows, cols);
  for(ni=0; ni<rows; ni++)
  {
     for(nj=0; nj<cols; nj++)
     {
        arr[ni][nj] =  element_value;
        printf("a[%d][%d] = %g\t", ni, nj, arr[ni][nj]; 
     }
     printf ("\n");
  }
}

// function to multiply 2-D matrices A, B and C using Blocked kij
void multipy_ABC_coded_temporal(double **arr_a, double **arr_b, double **arr_c, int rows, int cols, int block)
{
   int bi, bj, bk, ni, nj, nk;

   printf("  Blocked kij for C[bi+ni][bj+nj] ... \n");
   printf("     C[%d+%d][%d+%d]", bi, ni, bj, nj)
   for (bk=0; bk<rows; bk+=block) 
   {
      for (bi=0; bi<cols; bi+=block)
      {
	      for (bj=0; bj<rows; bj+=block) 
         {
			   for (nk=0; nk<block; nk++)
			   {
				   for (ni=0; ni<block; ni++) 
				   {
					   printf("   ");  
					   double x = arr_a[bi+ni][bk+nk];
					   for (nj=0; nj<block; nj++)
					   {
						   arr_c[bi+ni][bj+nj] += x * arr_b[bk+nk][bj+nj];
                     printf("     ", arr_c[bi+ni][bj+nj] );
					   }
                  printf("\n");
               }
            }
         }
      }
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

// start of program
   printf("\n");
   printf("\n......................................................................................................... \n");
   printf("\nSTART :\n");
 
   printf("INPUT : \n   Blocked kij (here nx = ny, i.e. |row| = |col|) ...\n");

// Get n value for row x col matrix, both equal to n
   printf("   Please enter matrix dimension --> n : ");
   int nx = get_value(MAXNVALUE);
   int ny = nx;

// Get block size
   printf("   Please enter sub-matrix dimension of n --> b : ");
   int nb = get_value(MAXNVALUE);
   while ( (!((nx % nb)==0)) && (!(nb<nx)) )
   {
      printf("   Error : (i) b not < n and (n mod b) not = 0. Please re-enter b : ");
      nb = get_value(MAXNVALUE);
   }

// Allocate memory for matrices A, B & C
   printf("CREATE MATRICES : \n");
   double **A = create_2Dmatrix(nx, ny);
   double **B = create_2Dmatrix(nx, ny);
   double **C = create_2Dmatrix(nx, ny);

// Initialize matrices A & B & C
   printf("INITIALIZE MATRICES : \n");
   init_2Dmatrix_random(A, nx, ny);
   init_2Dmatrix_random(B, nx, ny);
   init_2Dmatrix_specific(C, nx, ny, 0.0);

// Calculate program execution time taken : start
   gettimeofday(&tv1, &tz);

// Execute multiplication of blocked kij 2-D matrix NxN
   printf("RESULTS : Improving Spatial Locality - use a temporary variable \n");
   multipy_ABC_coded_spatial(A, B, C, nx, ny, nb);

// program execution time taken : end
   gettimeofday(&tv2, &tz);
   double elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

// output results, continued
   printf("   Time taken to complete : %g seconds.\n", elapsed);

// Calculate program execution time taken : start
   gettimeofday(&tv1, &tz);

// Execute multiplication of blocked kij 2-D matrix NxN
   printf("RESULTS : Improving Temporal Locality - blocked matrix multiplication \n");
   multipy_ABC_coded_temporal(A, B, C, nx, ny, nb);

// program execution time taken : end
  gettimeofday(&tv2, &tz);
  double elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

// output results, continued
  printf("   Time taken to complete : %g seconds.\n", elapsed);
	
// De-allocate memory for matrices A, B & C
  printf("CLEAN-UP : \n");
   deallocate_2Dmatrix_memory(A, nx);
   deallocate_2Dmatrix_memory(B, nx);
   deallocate_2Dmatrix_memory(C, nx);

//end of program
  printf("\nEND : \n");
  printf("\n......................................................................................................... \n");
  printf("\n");

  return 0;
}

