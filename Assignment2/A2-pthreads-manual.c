/* 
* *******************************************************************************

    Paula Dwan
    13208660 : COMP40730 : Assignment 4

    MPI program computing norm of product of two nxn matrices
    1/     Granularity of the product : 
           (ii)      One step algorithim. No intermediate resulting matrix
    2/     Partitioning Scheme :
           (A)      Left matrix is horizonitally partitioned
    3/     Matrix norm to be computed :
           (B)      Maximum absolute row sum norm (aka infinity-norm)
                     calculate max value of sum of each row (or product of same) 
                     ||A|| infintity = max SUM | Aij |
                     where   max     :==> 0 <= i <=n
                                  SUM   :==>  j =0 to n-1      

    The program A2-pthreads-solo.c calculates |C| and the infinity norm of |C| using :
    (i)     Straight-forward IJK 
    (ii)    pThreads 
    
    The program A2-pthreads-manual.c computes the time taked for manual evaulation of |C| and 
    the infinity norm of |C| using :
    (i)     Straight-forward IJK 
    (ii)    DGEMM computations.                                                          
                                  
* ********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>
#include <string.h>

void usage () 
{
    fprintf(stdout,"\nUSAGE :\t<program name> [<-r>|<-i>] [N] <matrix contents file>.txt <timing file>.dat \n");
    fprintf(stdout,"\nTO :\tCalculate |C| = |A| x |B| manually for pThreads comparison and also calculate infinity norm of |C|. \n");    
    fprintf(stdout,"\nWHERE :\t1. \t<-r> \tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n");
    fprintf(stdout," \t \t<-i> \tinitialize |A| & |B| _incrementally_ with <column> value and |C| with '0' \n");
    fprintf(stdout," \t2. \t[N] \tmax size of each matrix, if invalid defaults to 1,000 \n");
    fprintf(stdout," \t3. \t<matrix contents file>.txt\n \t \tname of .txt file to store values of matrices |A| |B| & |C| \n");
    fprintf(stdout," \t4. \t<timing .dat file> .dat \n \t \tname of .dat file to contain time to complete for each iteration \n \n");
    fprintf(stdout,"\tMANUAL \tStraight-forward IJK & DGEMM computations only \n");
    fprintf(stdout,"\tSOLO \tStraight-forward IJK & pThreads computations only \n \n");
    exit(0);
}

int validate_if_file_exists(char * fn)
{
    FILE *fp = fopen(fn, "r") ; 
    if (fp !=  NULL)
    {
        fclose(fp) ; 
        return 1 ;       // file exists
    }
    return 0 ;           // files does not exist
}

double* allocate_memory_matrix (int rows, int cols)
{
    double* l_matrix;  
    l_matrix = (double *) malloc(rows*cols*sizeof (double));  
    if(!(l_matrix))
    {
        fprintf(stderr,"ERROR :\texiting - memory allocation failed for matrix.\n");
        exit(1);
    }
    return l_matrix;  
} 

void deallocate_matrix_memory(double *l_matrix)
{
    free(l_matrix); 
}

void print_matrix(double *l_matrix, int rows, int cols, FILE *fp)
{
    int ni;
    int count =1;
    fprintf(fp, "\t\t");
    for (ni=0; ni<(rows*cols); ni++)
    {
        fprintf(fp, "%g \t", l_matrix[ni]); 
        if (count == rows) 
        {
            fprintf(fp, "\n\t\t");
            count =1;
        } else 
        {
            count ++;    
        }
    }
    fprintf(fp, "\n");    
}

void initialize_matrix_random(double *l_matrix, int rows, int cols)
{
    int ni, mod_rand=10;
    for (ni=0; ni<(rows*cols); ni++)
    {
        l_matrix[ni]= (double) ((rand() % mod_rand) + 1);
    }
}

void initialize_matrix_increment(double *l_matrix, int rows, int cols)
{
    int ni;    
    for (ni=0; ni<(rows*cols); ni++)
    {
        l_matrix[ni]= ((ni % rows) + 1);
    }
}

void initialize_matrix_zero(double *l_matrix, int rows, int cols)
{
    int ni;
    for (ni=0; ni<(rows*cols); ni++)
    {
        l_matrix[ni]=  0.0;
    }
}

void initialize_matrix_file_contents (FILE *fp) 
{
    fprintf(fp, "# --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA2-pthreads-manual \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "# --------------------------------------------------------------------------------------------------  \n ");
}

void initialize_data_file_contents (FILE *fp) 
{
    fprintf(fp, "# --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA2-pthreads-manual \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# |Matrix| \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \n# \n");
}

//  manual calculation
double multipy_abc_manual(double *matrix_A, double *matrix_B, double *matrix_C, int rows, int cols)
{
    int ni, nj, nk ; 
    
    for (ni=0 ; ni<rows ; ni++)
    {
        for (nj=0 ; nj<cols ; nj++)
        {
            double sum = 0.0 ; 
            for (nk=0 ; nk<rows ; nk++)
            {
                sum+= (matrix_A[(ni*rows)+nk]) * (matrix_B[(nk*rows)+nj]) ; 
            }
            matrix_C[(ni*rows)+nj] = sum ; 
        }
    }
double row_norm =0.0, inf_norm =0.0 ; 
    for (ni=0 ; ni<rows ; ni++)
    {
        row_norm =0.0 ; 
        for (nj=0 ; nj<rows ; nj++)
        {
            row_norm += matrix_C[(ni*rows) +nj] ; 
        }
        inf_norm = (inf_norm < row_norm) ? row_norm : inf_norm ; 
    }
    return inf_norm ; 
}

//  dgemm calculation 
double multipy_abc_dgemm(double *matrix_A, double *matrix_B, double *matrix_C, int rows, int cols)
{
    double row_norm =0.0, inf_norm =0.0 ; 
    int ALPHA=1.0 ; 
    int BETA=0.0 ; 
    int ni, nj ; 
//  m, n, k : local integers indicating the size of the matrices for rows x columns :: A :  m x  k, B :  k x n, C:  m x n
//  Here, m = n = k = rows = columns = <nx> = <ny> as supplied
    int lm = rows, ln = rows ; 
//  la_offset, lb_offset, lc_offset : Leading dimension of matrix A, B or C respectively, or the number of elements 
//  between successive rows for row-major storage or columns for column-major storage. 
    int la_offset = rows, lb_offset = cols, lc_offset = rows ; 
   
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lm, ln, ln, ALPHA, matrix_A, la_offset, matrix_B, lb_offset, BETA, matrix_C, lc_offset) ;   
    for (ni=0 ; ni<rows ; ni++)
    {
        row_norm =0.0 ; 
        for (nj=0 ; nj<rows ; nj++)
        {
            row_norm += matrix_C[(ni*rows) +nj] ; 
        }
        inf_norm = (inf_norm < row_norm) ? row_norm : inf_norm ; 
    }
    return inf_norm ; 
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main (int argc, char *argv[])
{

//  define variables
    struct timeval tv1, tv2;
    struct timezone tz;
    int increment_or_random = 0;    // increment enabled by default
    int max_num_args = 5;
    char filename_matrix[60];
    char filename_timing[60];
    int MAXN = 1000;

//  CLI PARAMETERS :: validate and initialize
    if (argc != max_num_args) 
    {
        fprintf(stderr, "\nERROR:\t<number of arguments> %d : is invalid, less than <default> %d\n", argc, max_num_args);      
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
    if (nx > MAXN)
    {
        fprintf(stderr, "\nWARNING:\tMatrix size entered <nx> [%d]  too large, now set to [%d]. \n", nx, MAXN); 
        nx = MAXN;
    }    
    int ny = nx;        
//  matrix file name .txt
    strncpy(filename_matrix, argv[3], 59);
    filename_matrix[60] = '\0';
    FILE *fp_matrix;
    int file_matrix_exists = validate_if_file_exists(filename_matrix);
    if (file_matrix_exists == 0) 
    {   
        fp_matrix= fopen(filename_matrix, "wa");
        initialize_matrix_file_contents(fp_matrix); 
    } else 
    {
        fp_matrix = fopen(filename_matrix, "a");
    } 
//  data file name .dat
    strncpy(filename_timing, argv[4], 59);
    filename_timing[60] = '\0';    
    FILE *fp_timing; 
    int file_timing_exists = validate_if_file_exists(filename_timing);    
    if (file_timing_exists == 0) 
    {   
        fp_timing= fopen(filename_timing, "wa");
        initialize_data_file_contents(fp_timing); 
    } else 
    {
        fp_timing = fopen(filename_timing, "a");
    }
    fprintf(fp_matrix, "# \n# RUNNING :\t%s %.2s %d \n", argv[0], init_type, nx);
    fprintf(stdout, "\n# RUNNING :\t%s %.2s %d \n", argv[0], init_type, nx);

//  CREATE & INITIALIZE :: matrices A & B & C and output results to matrix file for reference
    fprintf(stdout, "# ALLOCATE :\tmatrices |A|, |B| and |C| ... \n") ; 
    fprintf(fp_matrix, "\n# ALLOCATE :\tmatrices |A|, |B| and |C| ... \n") ; 
    double *A = allocate_memory_matrix(nx, ny);
    double *B = allocate_memory_matrix(nx, ny);
    double *C = allocate_memory_matrix(nx, ny);
    fprintf(stdout,"# INITIALIZE :\t|A| & |B| ... \n");
    fprintf(fp_matrix, "# INITIALIZE :\t|A| & |B| ... \n");
    if (increment_or_random == 1) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |A|  using random mod 10 ... \n", nx, ny);
        initialize_matrix_random(A, nx, ny);
        print_matrix(A, nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |B| using random mod 10 ... \n", nx, ny);
        initialize_matrix_random(B, nx, ny);
        print_matrix(B, nx, ny, fp_matrix);
    } else if (increment_or_random == 0) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |A| using incremental <column> value + 1 ... \n", nx, ny);
        initialize_matrix_increment(A, nx, ny);
        print_matrix(A, nx, ny, fp_matrix);            
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |B| using incremental <column> value + 1 ... \n", nx, ny);
        initialize_matrix_increment(B, nx, ny);
        print_matrix(B, nx, ny, fp_matrix);            
    }

//  MANUAL execution : calculate |C| and infinity norm for resulting |C|
    fprintf(stdout,"# INITIALIZE :\t<%d> x <%d> matrix |C| for Straight-forward IJK manual computation ... \n", nx, ny);
    fprintf(fp_matrix, "# INITIALIZE :\t<%d> x <%d> matrix |C| for Straight-forward IJK manual computation ... \n", nx, ny) ; 
    initialize_matrix_zero(C, nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double manual_norm = multipy_abc_manual(A, B, C, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double manual_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6 ; 
    fprintf(fp_matrix, "# RESULTS :\tmanual Straight-forward IJK calculation ... \n") ; 
    fprintf(fp_matrix, "# \t\tComputed Matrix [%d] x [%d] |C| ... \n", nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ;     
    fprintf(fp_matrix, "# \t\tMatrix |C| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", manual_elapsed, manual_norm) ; 
    fprintf(stdout, "# RESULTS :\tmanual Straight-forward IJK calculation ... \n") ; 
    fprintf(stdout, "# \t\tMatrix |C| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", manual_elapsed, manual_norm) ; 

//  CALCULATION :: |C| using DGEMM
    fprintf(stdout,"# INITIALIZE :\t|C| for BLAS/ATLAS computation ... \n");
    fprintf(fp_matrix, "# INITIALIZE :\t<%d> x <%d> matrix |C| for BLAS/ATLAS computation ... \n", nx, ny) ; 
    initialize_matrix_zero(C, nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double dgemm_norm = multipy_abc_dgemm(A, B, C, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6 ; 
    fprintf(stdout, "# RESULTS :\tBLAS/ATLAS computation ...\n") ; 
    fprintf(stdout, "# \t\tMatrix |C| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", dgemm_elapsed, dgemm_norm) ;     
    fprintf(fp_matrix, "\n# RESULTS :\tBLAS/ATLAS computation ...\n") ; 
    fprintf(fp_matrix, "\n# \t\tComputed Matrix [%d] x [%d] |C| ... \n", nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ; 
    fprintf(fp_matrix, "# \t\tMatrix |C| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", dgemm_elapsed, dgemm_norm) ; 

//  OUTPUT :: results to stdout & .dat file : |Matrix| ||inf norm/manual || Time/manual || infinity norm || Time / dgemm
    fprintf(stdout,"# SUMMARY :\t|Matrix|  \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \n");
    fprintf(stdout,"# \t\t%d \t\t%f \t%.1f \t\t%f \t%.1f \n", nx, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm);
    fprintf(fp_timing, "%d \t%f \t%.1f \t%f \t%.1f \n", nx, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm);
   	
//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(A);
    deallocate_matrix_memory(B);
    deallocate_matrix_memory(C);    
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
