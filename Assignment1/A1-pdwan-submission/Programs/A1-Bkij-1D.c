/* 
********************************************************************************

    Paula Dwan
    13208660 : COMP40730 : Assignment 1

    Blocked KIJ IMPLEMENTATION

    Write C programs implementing the following three algorithms of multiplication of two n×n 
    dense matrices:
    1)    Straightforward non-Blocked ijk algorithm.
    2)    Blocked ijk algorithm using square b×b blocks.
    3)    Blocked kij algorithm using square b×b blocks.

    Experiment with the programs and build/plot:
    1)  The dependence of the execution time of each program on the matrix size n and the 
          block size b .
          b.  BLAS calls
    2)  The speedup of the Blocked algorithms over the non-Blocked one as a function of the  
          matrix size and the block size.
          a.  manually written      
    3)  Compare the fastest program with the BLAS/ATLAS routine dgemm implementing the 
          same operation. 
          Explain the results.
          c.  ATLAS      
   
*********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>
#include <string.h>

void usage () 
{
    fprintf(stdout,"\nUSAGE : \t<program name> [<-r>|<-i>] [N] <matrix contents file>.txt <data timing file>.dat \n");
    fprintf(stdout,"\nTO : \t\tCalculate |C| = |A| x |B| using algorithm : Blocked KIJ using square bxb block \n");    
    fprintf(stdout,"\nWHERE :");
    fprintf(stdout,"\t1. \t<-r>\tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n");
    fprintf(stdout,"\t   \t<-i>\tinitialize |A| & |B| _incrementally_ with <column> value and |C| with '0' \n");
    fprintf(stdout,"\t2. \t[N] \tmax size of each matrix, if invalid (greater than maximum set), then set to [${MAXN}]. \n");
    fprintf(stdout,"\t3. \t[B] \tblock size appliable to all matrices, if invalid (i) greater than max permitted or (ii) [N] % [B] not equal zero.  \n\t\t\tSet to [${MAXB}].\n");
    fprintf(stdout,"\t4. \t<matrix contents file>.txt\n\t\tname of .txt file to store values of matrices |A| |B| & |C| \n");
    fprintf(stdout,"\t5. \t<data timingfile>.dat \n\t\tname of .dat file to contain time to complete for each iteration \n\n");
    exit(0);
}

double* allocate_memory_matrix (int rows, int cols)
{
    double* l_matrix;  
    l_matrix = (double *) malloc(rows*cols*sizeof (double));  
    if(!(l_matrix))
    {
        fprintf(stderr,"ERROR : \texiting - memory allocation failed for matrix.\n");
        exit(1);
    }
    return l_matrix;  
} 

void deallocate_matrix_memory(double *matrix1d)
{
    free(matrix1d); 
}

void print_matrix(double *matrix1d, int rows, int cols, FILE *fp)
{
    int ni;
    int count =1;
    for (ni=0; ni<(rows*cols); ni++)
    {
        fprintf(fp, "%g\t", matrix1d[ni]); 
        if (count == rows) 
        {
            fprintf(fp, "\n");
            count =1;
        } else 
        {
            count ++;    
        }
    }
}

void init_matrix_random(double *matrix1d, int rows, int cols)
{
    int ni, mod_rand=10;
    for (ni=0; ni<(rows*cols); ni++)
    {
        matrix1d[ni]= (double) ((rand()%mod_rand) + 1);
    }
}

void init_matrix_increment(double *matrix1d, int rows, int cols)
{
    int ni;    
    for (ni=0; ni<(rows*cols); ni++)
    {
        matrix1d[ni]= ((ni % rows) + 1);
    }
}

void init_matrix_zero(double *matrix1d, int rows, int cols)
{
    int ni;
    for (ni=0; ni<(rows*cols); ni++)
    {
        matrix1d[ni]=  0.0;
    }
}

int validate_if_file_exists(char * fn)
{
    FILE *fp = fopen(fn, "r") ;
    if ( fp !=  NULL )
    {
        fclose(fp);
        return 1;       // file exists
    }
    return 0;           // files does not exist
}

void init_matrix_file_contents(FILE *fp) 
{
    fprintf(fp, "#  -------------------------------------------------------------------------------------------------... \n# \n");
    fprintf(fp, "# Program :\tA1-Bkij-1D\n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n");
}

void init_timing_file_contents(FILE *fp) 
{
    fprintf(fp, "#  -------------------------------------------------------------------------------------------------... \n# \n");
    fprintf(fp, "# Program :\tA1-Bkij-1D\n# where :\t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "#  -------------------------------------------------------------------------------------------------... \n");
    fprintf(fp, "# |Matrix| \t|Block| \tTime/straight-forward \tTime/Blocked \tTime/dgemm \n# \n");
}

void multipy_ABC_cblas(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols, int block, double alpha, double beta)
{
//                  m, n, k :    local integers indicating the size of the matrices for rows x columns :: A :  m x  k, B :  k x n, C:  m x n. Here, m = n = k = block
//                  la_offset, lb_offset, lc_offset : Leading dimension of matrix A, B or C respectively, or the number of elements between 
//                 successive rows for row-major storage or columns for column-major storage. 
    int lm = block, ln = block; 
    int la_offset = rows, lb_offset = cols, lc_offset = rows;
    int i, j, k;
    for (k = 0; k < lc_offset ; k += block)
    {
        for (i = 0; i < la_offset; i += block)
        {
            for (j = 0; j < lb_offset; j += block)
            {
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lm, ln, ln, alpha, matrix_a+i*la_offset+k, la_offset, matrix_b+j+k*lb_offset, lb_offset, beta, matrix_c+i*lc_offset+j, lc_offset);
            }
        }
    }
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main ( int argc, char *argv[] )
{

//  define variables
    struct timeval tv1, tv2;
    struct timezone tz;
    const double ALPHA = 1.0;
    const double BETA = 1.0;
    int increment_or_random = 5;
    int max_num_args = 6;
    char filename_matrix[50];
    char filename_timing[50];
    int MAXN = 1000;
    int MAXB=100;
    int ni, nj, nk, bi, bj, bk;
    double simple_elapsed=0.0, complex_elapsed=0.0, dgemm_elapsed=0.0;
          
//  CLI PARAMETERS :: validate and initialize
    if ( argc != max_num_args ) 
    {
        fprintf(stderr, "ERROR : \t<number of arguments> : %d, is invalid, less than <default> : %d. \n", argc, max_num_args);        
        usage();
    }
//  random or increment initialization of matrices |A| and |B|
    char init_type[3];
    strncpy(init_type, argv[1], 2);
    init_type[3] = '\0';
    if (strncmp(init_type,"-i", 2) == 0) 
    {   
        increment_or_random = 0; 
    } else if (strncmp(init_type,"-r", 2) == 0) 
    {   
        increment_or_random = 1; 
    } else
    { 
        fprintf(stderr, "ERROR : \t'invalid entry : %s for '-i' or '-r'. \n", init_type);        
        usage();
    }
//  matrix size
    int nx = atoi(argv[2]);                                     
    if ( nx > MAXN )
    {
        fprintf(stderr, "WARNING : \tMatrix size entered <nx> %d too large, now set to %d. \n", nx, MAXN);
        nx = MAXN;
    }    
    int ny = nx;  
//  block size 
    int nb = atoi(argv[3]);                                     
    if ( nb > MAXB )
    {
        fprintf(stderr, "WARNING : \tMatrix size entered <nb> %d  too large, now set to %d. \n", nb, MAXB);
        nb = MAXB;        
    } else if (nb > nx)
    {
        fprintf(stderr, "ERROR : \t<nb> %d : Block size must be less <nx> %d : matrix size. \n", nb, nx);
        usage();
    } else if ( (nx % nb) != 0 )
    {
        fprintf(stderr, "WARNING : \tBlock size <nb> %d is not an even multiple of Matrix size <nx> %d. Both set to defaults of <MAXN> %d and <MAXB> %d.\n", nb, nx, MAXN, MAXB);
        nx = MAXN;
        nb = MAXB;
    }    
//  matrix file name .txt
    strncpy(filename_matrix, argv[4], 49);
    filename_matrix[50] = '\0';
    FILE *fp_matrix;
    int file_matrix_exists = validate_if_file_exists(filename_matrix);
    if ( file_matrix_exists == 0 ) 
    {   
        fp_matrix= fopen(filename_matrix, "wa" );
        init_matrix_file_contents(fp_matrix); 
    } else 
    {
        fp_matrix = fopen(filename_matrix, "a+" );
    } 
//  data file name .dat
    strncpy(filename_timing, argv[5], 49);
    filename_timing[50] = '\0';    
    FILE *fp_timing; 
    int file_timing_exists = validate_if_file_exists(filename_timing);
    if ( file_timing_exists == 0 ) 
    {   
        fp_timing= fopen(filename_timing, "wa" );
        init_timing_file_contents(fp_timing); 
    } else 
    {
        fp_timing = fopen(filename_timing, "a+" );
    } 
    fprintf(fp_matrix, "# \n# RUNNING : \t%s %.2s %d %d \n", argv[0], init_type, nx, nb);
    fprintf(stdout, "\n# RUNNING : \t%s %.2s %d %d \n", argv[0], init_type, nx, nb);

//  CREATE & INITIALIZE :: matrices A & B & C and output results to matrix file for reference
    fprintf(stdout,"# CREATE MATRICES ... \n");
    double *A =  allocate_memory_matrix(nx, ny);
    double *B =  allocate_memory_matrix(nx, ny);
    double *C =  allocate_memory_matrix(nx, ny);

    fprintf(stdout,"# INITIALIZE MATRICES ... \n");
    if (increment_or_random == 1) 
    {
        init_matrix_random(A, nx, ny);
        init_matrix_random(B, nx, ny);
    } else if (increment_or_random == 0) 
    {
        init_matrix_increment(A, nx, ny);
        init_matrix_increment(B, nx, ny);
    }
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |A| ... \n", nx, ny);
    print_matrix(A, nx, ny, fp_matrix);
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |B| ... \n", nx, ny);
    print_matrix(B, nx, ny, fp_matrix);

//  CALCULATION :: |C| using Straight-forward KIJ manual computation
    fprintf(stdout,"# RESULTS : manual straightforward KIJ manual computation ... \n");
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C| for manual Straight-forward KIJ computation ... \n", nx, ny);
    init_matrix_zero(C, nx, ny);
    print_matrix(C, nx, ny, fp_matrix);       
    gettimeofday(&tv1, &tz);
    for (nk=0; nk<nx; nk++)
    {
        for (ni=0; ni<ny; ni++)
        {
            for (nj=0; nj<nx; nj++)
            {
                C[(nk*nx)+ni] += (A[(nk*nx)+nj]) * (B[(nj*nx)+ni]);
            }
        }
    }
    gettimeofday(&tv2, &tz);
    simple_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# |C| : <%d> x <%d> matrix computed values for manual Blocked KIJ computation ... \n", nx, ny);
    print_matrix(C, nx, ny, fp_matrix);
    fprintf(fp_matrix, "# |C| : matrix calculated in %f seconds ... \n", simple_elapsed);


//  CALCULATION :: |C| using Blocked KIJ manual computation
    fprintf(stdout,"# RESULTS : manual Blocked KIJ computation ... \n");
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C| for manual Blocked KIJ computation .. \n", nx, ny);
    init_matrix_zero(C, nx, ny);
    print_matrix(C, nx, ny, fp_matrix);       
    gettimeofday(&tv1, &tz);
    for(bk = 0; bk < nx; bk += nb) 
    {
      for(bi = 0; bi < nx; bi += nb) 
      {
         for(bj = 0; bj < nx; bj += nb) 
         {
            for(nk = 0; nk < nb; nk++) 
            {
               for(ni = 0; ni < nb; ni++) 
               {
                  for(nj = 0; nj < nb; nj++) 
                  {
                        C[ (bk + nk) * nx + bi + ni] += A[ (bk + nk) * nx + bi + nj] * B[ (bi + nj) * nx + bi + ni];
                        // fprintf(stdout, "DEBUG : |C|  += |A| * |B| : C[%d] +=A[%d] * C[%d] : %g += %g * %g \n",(bk+nk)*nx +bi+ni, (bk+nk)*nx +bi+nj,  (bi+nj) *nx +bi+ni, C[ (bk+nk)*nx +bi+ni],  A[ (bk+nk)*nx +bi+nj], B[ (bi+nj) *nx +bi+ni] );
                  }
               }
            }
         }
      }
    }
    gettimeofday(&tv2, &tz);
    complex_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# |C| : <%d> x <%d> matrix computed values for CBLAS/ATLAS Blocked KIJ computation ... \n", nx, ny);
    fprintf(fp_matrix, "# |C| : matrix calculated in %f seconds ... \n", complex_elapsed);
    print_matrix(C, nx, ny, fp_matrix);

//  CALCULATION :: |C| using DGEMM
    fprintf(stdout,"# RESULTS : BLAS/ATLAS KIJ computation ... \n");
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C| for CBLAS/ATLAS Blocked KIJ computation ... \n", nx, ny);
    init_matrix_zero(C, nx, ny);
    print_matrix(C, nx, ny, fp_matrix);    
    gettimeofday(&tv1, &tz);
    multipy_ABC_cblas(A, B, C, nx, ny, nb, ALPHA, BETA);
    gettimeofday(&tv2, &tz);
    dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# |C| : <%d> x <%d> matrix computed values using CBLAS/ATLAS KIJ Blocked KIJ computation ... \n", nx, ny);
    fprintf(fp_matrix, "# |C| : calculated in %f seconds... \n",  dgemm_elapsed);
    print_matrix(C, nx, ny, fp_matrix);

//  OUTPUT :: results to stdout & .dat file : matrix size || Time/manual || Time / dgemm
    fprintf(stdout,"# \t\t|Matrix|   |Block|   Time/straight-forward   Time/Blocked   Time/dgemm\n");
    fprintf(stdout,"# Results: \t%d \t%d \t%f \t%f \t%f \n", nx, nb, simple_elapsed, complex_elapsed, dgemm_elapsed);
    fprintf(fp_timing, "%d \t%d \t%f \t%f \t%f \n", nx, nb, simple_elapsed, complex_elapsed, dgemm_elapsed);
   	
//  CLEANUP :: deallocate memory & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(A);
    deallocate_matrix_memory(B);
    deallocate_matrix_memory(C);    
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
