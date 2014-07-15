/* 
********************************************************************************
    STRAIGHT-FORWARD IJK IMPLEMENTATION

    Write C programs implementing the following three algorithms of multiplication of two n×n 
    dense matrices:
    1)    Straightforward non-blocked ijk algorithm.        <*****
    2)    Blocked ijk algorithm using square b×b blocks.
    3)    Blocked kij algorithm using square b×b blocks.

    Experiment with the programs and build/plot:
    1)  The dependence of the execution time of each program on the matrix size n and the 
          block size b .
          b.  BLAS calls
    2)  The speedup of the blocked algorithms over the non-blocked one as a function of the  
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
    fprintf(stdout,"\nTO : \t\tCalculate |C| = |A| x |B| using algorithm : Straightforward IJK \n");    
    fprintf(stdout,"\nWHERE :");
    fprintf(stdout,"\t1. \t<-r>\tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n");
    fprintf(stdout,"\t   \t<-i>\tinitialize |A| & |B| _incrementally_ with <column> value and |C| with '0' \n");
    fprintf(stdout,"\t2. \t[N] \tmax size of each matrix, if invalid (greater than max), then defaults to [${MAXN}].\n");
    fprintf(stdout,"\t3. \t<matrix contents file>.txt\n\t\tname of .txt file to store values of matrices |A| |B| & |C| \n");
    fprintf(stdout,"\t4. \t<data timingfile>.dat \n\t\tname of .dat file to contain time to complete for each iteration \n\n");
    exit(0);
}

double* allocate_memory_matrix (int rows, int cols)
{
    double* l_matrix;  
    l_matrix = (double *) malloc(rows*cols*sizeof (double));  
    if(!(l_matrix))
    {
        fprintf(stderr,"ERROR : \texiting - memory allocation failed for matrix.\n");
        usage();
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
        return 1;       // exists
    }
    return 0;          // does not exist
}

void init_matrix_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA1-Sijk-1D\n# where :\t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n# \n");
}

void init_timing_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA1-Sijk-1D\n# where :\t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "# Data values from each run using different matrix and block size. \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# |Matrix|   Time/simple  Time/complex  Time/dgemm \n# \n");
}

void multipy_ABC_cblas(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols, double alpha, double beta)
{
// m, n, k :    local integers indicating the size of the matrices for rows x columns :: A :  m x  k, B :  k x n, C:  m x n
//                 Here, m = n = k = rows = columns = <nx> = <ny> as supplied
    int lm = rows, ln = rows;
// la_offset, lb_offset, lc_offset :
//                 Leading dimension of matrix A, B or C respectively, or the number of elements between 
//                 successive rows for row-major storage or columns for column-major storage. 
    int la_offset = rows, lb_offset = cols, lc_offset = rows;
    cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, lm, ln, ln, alpha, matrix_a, la_offset, matrix_b, lb_offset, beta, matrix_c, lc_offset);   
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main ( int argc, char *argv[] )
{

//  define variables
    struct timeval tv1, tv2;
    struct timezone tz;
    const double ALPHA = 1.0;
    const double BETA = 0.0;
    int increment_or_random = 5;
    int max_num_args = 5;
    char filename_matrix[50];
    char filename_timing[50];
    int MAXN = 1000;
    int ni, nj, nk;
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
        fprintf(stderr, "\nERROR : \t'invalid entry : %s for '-i' or '-r'", init_type);        
        usage();
    }
    //  matrix size
    int nx = atoi(argv[2]);                                     
    if ( nx > MAXN )
    {
        fprintf(stderr, "WARNING : \tMatrix size entered <nx> %d too large, now set to %d \n", nx, MAXN);
        nx = MAXN;
    }    
    int ny = nx;  
//  matrix file name .txt
    strncpy(filename_matrix, argv[3], 49);
    filename_matrix[50] = '\0';
    FILE *fp_matrix;
    int file_matrix_exists = validate_if_file_exists(filename_matrix);
    if ( file_matrix_exists == 0 ) 
    {   
        fp_matrix= fopen(filename_matrix, "wa" );
        init_matrix_file_contents(fp_matrix); 
    } else 
    {
        fp_matrix = fopen(filename_matrix, "a" );
    } 
//  data file name .dat
    strncpy(filename_timing, argv[4], 49);
    filename_timing[50] = '\0';    
    FILE *fp_timing; 
    int file_timing_exists = validate_if_file_exists(filename_timing);
    if ( file_timing_exists == 0 ) 
    {   
        fp_timing= fopen(filename_timing, "wa" );
        init_timing_file_contents(fp_timing); 
    } else 
    {
        fp_timing = fopen(filename_timing, "a" );
    } 
    fprintf(fp_matrix, "# \n# RUNNING : \t%s %.2s %d \n\n", argv[0], init_type, nx);
    fprintf(stdout, "# \n# RUNNING : \t%s %.2s %d \n\n", argv[0], init_type, nx);

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
    init_matrix_zero(C, nx, ny);

    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |A| ... \n", nx, ny);
    print_matrix(A, nx, ny, fp_matrix);
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |B| ... \n", nx, ny);
    print_matrix(B, nx, ny, fp_matrix);
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C| ... \n", nx, ny);
    print_matrix(C, nx, ny, fp_matrix);

//  CALCULATION :: |C| using manual simple
    fprintf(stdout,"# RESULTS : simple manual calculation ... \n");
    gettimeofday(&tv1, &tz);
    for (ni=0; ni<nx; ni++)
    {
        for (nj=0; nj<ny; nj++)
        {
            for (nk=0; nk<nx; nk++)
            {
                C[(ni*nx)+nj] += (A[(ni*nx)+nk]) * (B[(nk*nx)+nj]);
            }
        }
    }
    gettimeofday(&tv2, &tz);
    simple_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# |C| : <%d> x <%d> matrix computed values : MANUAL simple ... \n", nx, ny);
    print_matrix(C, nx, ny, fp_matrix);
    fprintf(fp_matrix, "# |C| : matrix calculated in %f seconds ... \n", simple_elapsed);

//  CALCULATION :: |C| using manual complex
    fprintf(stdout,"# RESULTS : complex manual calculation ... \n");
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C|, redone for MANUAL complex .. \n", nx, ny);
    init_matrix_zero(C, nx, ny);
    print_matrix(C, nx, ny, fp_matrix);       
    gettimeofday(&tv1, &tz);
    for (ni=0; ni<nx; ni++)
    {
        for (nj=0; nj<ny; nj++)
        {
            double sum = 0.0;
            for (nk=0; nk<nx; nk++)
            {
                sum+= (A[(ni*nx)+nk]) * (B[(nk*nx)+nj]);
            }
            C[(ni*nx)+nj] = sum;
        }
    }
    gettimeofday(&tv2, &tz);
    complex_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# |C| : <%d> x <%d> matrix computed values : MANUAL complex ... \n", nx, ny);
    print_matrix(C, nx, ny, fp_matrix);
    fprintf(fp_matrix, "# |C| : matrix calculated in %f seconds ... \n", complex_elapsed);

//  CALCULATION :: |C| using DGEMM
    fprintf(stdout,"# RESULTS : BLAS/ATLAS calculation -\n");
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C|, redone for CBLAS/ATLAS ... \n", nx, ny);
    init_matrix_zero(C, nx, ny);
    print_matrix(C, nx, ny, fp_matrix);    
    gettimeofday(&tv1, &tz);
    multipy_ABC_cblas(A, B, C, nx, ny, ALPHA, BETA);
    gettimeofday(&tv2, &tz);
    dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "\n# |C| : <%d> x <%d> matrix computed values using CBLAS/ATLAS ... \n", nx, ny);
    print_matrix(C, nx, ny, fp_matrix);
    fprintf(fp_matrix, "# |C| : calculated in %f seconds... \n",  dgemm_elapsed);

//  OUTPUT :: results to stdout & .dat file : matrix size || Time/simple || Time/complex || Time / dgemm
    fprintf(stdout,"# \t\t|Matrix|  Time/simple  Time/complex  Time/dgemm\n");
    fprintf(stdout,"# Results: \t%d \t%f \t%f \t%f \n", nx, simple_elapsed, complex_elapsed, dgemm_elapsed);
    fprintf(fp_timing, "%d \t%f \t%f \t%f \n", nx, simple_elapsed, complex_elapsed, dgemm_elapsed);
   	
//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(A);
    deallocate_matrix_memory(B);
    deallocate_matrix_memory(C);    
    fclose(fp_matrix);
    fclose(fp_timing);
    
    return 0;
}
