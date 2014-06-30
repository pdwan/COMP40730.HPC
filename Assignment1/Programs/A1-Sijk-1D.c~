/* ********************************************************************************
*
*   Paula Dwan
*   13208660 : COMP40730 : Assignment 1 : Part I
*
*   STRAIGHT-FORWARD IJK IMPLEMENTATION
* 
*   Write C programs implementing the following three algorithms of multiplication of two n×n 
*   dense matrices:
*   1)    Straightforward non-blocked ijk algorithm.        <*****
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
#include <string.h>


void usage () {
    printf( "USAGE : \t<program name> <-i | -r> <N> <matrix contents .txt file> <timing .dat file>");
    printf( "where ");
    printf("-r| \tinitialize A| & |B| with random numbers and |C| with '0' ");
    printf("\t-t \tinitialize |A| & |B| incrementally with <row> value and |C| with '0' ");
    printf("\t-t \tmax size of each matrix, defaults to 1000 if invalid - verified in calling script");
    printf( "\t<matrix contents .txt file>\n\t\tname of .txt file to contain the matrix contents");
    printf( "\t<timing .dat file>\n\t\tname of .dat file to contain time to complete for each iteration");
    exit(1);
}

void add_comments_to_files (FILE *fp) 
{
    fprintf(fp, "# ......................................................................................................... \n");
    fprintf(fp, "# Computation of Blocked IJK using square bxb block for matrices |A| & |B| storing results in |C|") ;
    fprintf(fp, "#\n# Matrix Size\t Block Size\t Time/manual\t Time/manual\t Time/dgenn \n #\t\tSimple\tComplex");	
    fprintf(fp, "# ......................................................................................................... \n");
}

double* allocate_memory_1d_matrix (int rows, int cols)
{
    double* l_matrix;  
    l_matrix = (double *) malloc(rows*cols*sizeof (double));  
    if(!(l_matrix))
    {
        printf("ERROR : \texiting - memory allocation failed for matrix.\n");
        exit(1);
    }
    return l_matrix;  
} 

void deallocate_1d_matrix_memory(double *matrix1d)
{
    free(matrix1d); 
}

void init_1d_matrix_random(double *matrix1d, int rows, int cols, FILE *fp)
{
    int ni, mod_rand=10;
    int count =1;
    fprintf(fp, "\n|A| or |B| : \tMatrix initialization using (mod %d)+1 on random number for matrix of size <%d> x <%d> ...\n", mod_rand, rows, cols);
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
    fprintf(fp, "\n|A| or |B| : \tMatrix initialization with row value for matrix of size <%d> x <%d> ...\n", rows, cols);
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
    fprintf(fp, "\n|C| : \tMatrix initialization with value of zero for matrix of size <%d> x <%d> ...\n", rows, cols);
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

void multipy_ABC_simple(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols, FILE *fp)
{
   int ni, nj, nk;
   int count = 1;
   for (ni=0; ni<rows; ni++)
   {
      for (nj=0; nj<cols; nj++)
      {
         for (nk=0; nk<rows; nk++)
         {
            matrix_c[ni+(nj*cols)] += (matrix_a[(ni*rows)+nk]) * (matrix_b[(nk*rows)+nj]);
         }
         fprintf(fp, "%g\t",matrix_c[ni+(nj*cols)]); 
         if (count == rows) {
            fprintf(fp, "\n");
            count =1;
        }else {
        count++;
        }
      }
   }
}

void multipy_ABC_complex(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols, FILE *fp)
{
   int ni, nj, nk;
   int count = 1;
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
         fprintf(fp, "%g\t",matrix_c[ni+(nj*cols)]); 
         if (count == rows) {
            fprintf(fp, "\n");
            count =1;
        }else {
        count++;
        }
      }
   }
}

void multipy_ABC_cblas(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols, double alpha, double beta, FILE *fp)
{
   int ni;
   int count =1;
// m, n, k :    local integers indicating the size of the matrices for rows x columns :: A :  m x  k, B :  k x n, C:  m x n
//                 Here, m = n = k = rows = columns = <nx> = <ny> as supplied
   int lm = rows, ln = rows;
// la_offset, lb_offset, lc_offset :
//                 Leading dimension of matrix A, B or C respectively, or the number of elements between 
//                 successive rows for row-major storage or columns for column-major storage. 
   int la_offset = rows, lb_offset = cols, lc_offset = rows;

   cblas_dgemm( CblasRowMajor, CblasTrans, CblasNoTrans, lm, ln, ln, alpha, matrix_a, la_offset, matrix_c, lb_offset, beta, matrix_c, lc_offset);   
    for(ni=0; ni<(rows*cols); ni++)
    {
        fprintf(fp, "%g\t", matrix_c[ni]); 
        if (count == rows) {
            fprintf(fp, "\n");
            count =1;
        }else {
        count++;
        }
    }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main ( int argc, char *argv[] )
{
// define variables
   struct timeval tv1, tv2;
   struct timezone tz;
   const double ALPHA = 1.0;
   const double BETA = 1.0;
   int increment_or_random = 5;
   int max_num_args = 5;

//  USAGE: <program name> <-i | -r> <N> <matrix contents .txt file> <timing .dat file>

    if ( argc != max_num_args )
    {
        usage();
    }

//  initialise values using parameters provided
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
    char filename_matrix[50];
    strncpy(filename_matrix, argv[3], 49);
    filename_matrix[50] = '\0';
    FILE *fp_matrix = fopen( filename_matrix, "wa" );
    char filename_timing[50];
    strncpy(filename_timing, argv[4], 49);
    filename_timing[50] = '\0';
    FILE *fp_timing = fopen(filename_timing, "wa" );
            
    fprintf(fp_matrix, "RUNNING : \t%s %s %d %s %s \n", argv[0],cell_value_type,nx,filename_matrix,filename_timing );

//  add initial comments to matrix values and timing data files 
    add_comments_to_files(fp_matrix);			
    add_comments_to_files(fp_timing);			

//  Allocate memory for matrices A, B & C
    fprintf(fp_matrix, "CREATE MATRICES : \n");
    double *A =  allocate_memory_1d_matrix(nx, ny);
    double *B =  allocate_memory_1d_matrix(nx, ny);
    double *C =  allocate_memory_1d_matrix(nx, ny);

//  Initialize matrices A & B & C
    fprintf(fp_matrix, "INITIALIZE MATRICES : \n");
    if (increment_or_random == 0) {
        init_1d_matrix_random(A, nx, ny, fp_matrix);
        init_1d_matrix_random(B, nx, ny, fp_matrix);
    } else if (increment_or_random == 1) 
    {
        init_1d_matrix_increment(A, nx, ny, fp_matrix);
        init_1d_matrix_increment(B, nx, ny, fp_matrix);
    }
    init_C_1d_matrix(C, nx, ny, fp_matrix);

//  MATRIX : simple matrix calculation
    gettimeofday(&tv1, &tz);
    fprintf(fp_matrix, "# RESULTS : simple manual calculation - \n");
    multipy_ABC_simple(A, B, C, nx, ny, fp_matrix);
    gettimeofday(&tv2, &tz);
    double elapsed_simple = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

//  MATRIX : complex matrix calculation
    gettimeofday(&tv1, &tz);
    fprintf(fp_matrix, "# RESULTS : simple manual calculation - \n");
    multipy_ABC_complex(A, B, C, nx, ny, fp_matrix);
    gettimeofday(&tv2, &tz);
    double elapsed_complex = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

    gettimeofday(&tv1, &tz);
    fprintf(fp_matrix, "# RESULTS : BLAS/ATLAS calculation -\n");
    multipy_ABC_cblas(A, B, C, nx, ny, ALPHA,BETA, fp_matrix);
    gettimeofday(&tv2, &tz);
    double elapsed_cblas = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

//  output results to .dat file : matrix size || Time/manual 	|| Time / dgemm
    fprintf(fp_timing, "%d \t %g \t %g \t %g \n",  nx, elapsed_simple, elapsed_complex, elapsed_cblas);
   	
//  De-allocate memory for matrices A, B & C
    printf("CLEAN-UP : \n");
    deallocate_1d_matrix_memory(A);
    deallocate_1d_matrix_memory(B);
    deallocate_1d_matrix_memory(C);
    
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
