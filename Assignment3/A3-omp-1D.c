/* 
********************************************************************************
    Paula Dwan
    13208660 : COMP40730 : Assignment 3

    Open MP program computing norm of product of two nxn matrices
    1/     Granularity of the product : 
           (ii) One step algorithim. No intermediate resulting matrix
    2/     Partitioning Scheme :
           (a)  Left matrix is horizonitally partitioned
    3/     Matrix norm to be computed :
           (b)  Maximum absolute row sum norm (aka infinity-norm)
                calculate max value of sum of each row (or product of same) 
                    ||A|| infintity = max SUM | Aij |
                    where   max :==> 0 <= i <=n
                            SUM :==> j =0 to n-1      
*********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>
#include <string.h>
#include <omp.h>

void usage () 
{
    fprintf(stdout,"\nUSAGE : \t<program name> [<-r>|<-i>] [N] [T] <matrix contents file>.txt <timing file>.dat \n");
    fprintf(stdout,"\nTO : \t\tCalculate |C| = |A| x |B| using Open MP and also calculate infinity norm of |C| \n");    
    fprintf(stdout,"\nWHERE :");
    fprintf(stdout,"\t1. \t<-r>\tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n");
    fprintf(stdout,"\t   \t<-i>\tinitialize |A| & |B| _incrementally_ with <column> value and |C| with '0' \n");
    fprintf(stdout,"\t2. \t[N] \tmax size of each matrix, if invalid defaults to 1,000 \n");
    fprintf(stdout,"\t3. \t[T] \tnumber of threads (i) less than [N] and (ii) [N] mod [T] = 0 \n");
    fprintf(stdout,"\t4. \t<matrix contents file>.txt\n\t\tname of .txt file to store values of matrices |A| |B| & |C| \n");
    fprintf(stdout,"\t5. \t<timing .dat file> .dat \n\t\tname of timing data file to containing calculation time for each iteration \n\n");
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
    FILE *fp;
    if ( fp = fopen(fn, "r") )
    {
        fclose(fp);
        return 1;
    }
    return 0;
}

void init_matrix_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA3-omp-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n");
}

void init_timing_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA3-omp-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n");
    fprintf(fp, "# |Matrix| \t|Threads| \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \n# \n");
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
    int max_num_args = 6;
    char filename_matrix[50];
    char filename_timing[50];
    int MAXTHREADS = omp_get_num_threads();
    int MAXN = 1000;
    int ni, nj, nk, number_of_row_slices, omp_nt;
    double manual_row_norm=0.0, manual_inf_norm=0.0, dgemm_row_norm=0.0, dgemm_inf_norm=0.0;
    double manual_elapsed=0.0, dgemm_elapsed=0.0;
          
//  CLI PARAMETERS :: validate and initialize
    if ( argc != max_num_args ) 
    {
        fprintf(stderr, "\nERROR : \t<number of arguments> : %d, is invalid, less than <default> : %d. \n", argc, max_num_args);        
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
        fprintf(stderr, "\nWARNING : \tMatrix size entered <nx> %d  too large, now set to %d \n", nx, MAXN);
        nx = MAXN;
    }    
    int ny = nx;
//  number of threads      
    int nt = atoi(argv[3]);  
    if ( (nx % nt) != 0)
    {
        fprintf(stderr, "\nERROR : \t<nt> %d : number of threads must divide evenly into <nx> %d : matrix size \n", nt, nx);
        usage();
    }
    if (nt > nx)
    {
        fprintf(stderr, "\nERROR : \t<nt> %d : number of threads must be less <nx> %d : matrix size \n", nt, nx);
        usage();
    }
    omp_set_num_threads(nt);
    number_of_row_slices = nx / nt;
//  matrix file name .txt
    strncpy(filename_matrix, argv[4], 49);
    filename_matrix[50] = '\0';
    FILE *fp_matrix;
    int file_matrix_exists = validate_if_file_exists(filename_matrix);
    if ( file_matrix_exists == 0 ) 
    {   
        fp_matrix= fopen(filename_matrix, "w" );
        init_matrix_file_contents(fp_matrix); 
    } else 
    {
        fp_matrix = fopen(filename_matrix, "a " );
    } 
//  data file name .dat
    strncpy(filename_timing, argv[5], 49);
    filename_timing[50] = '\0';    
    FILE *fp_timing; 
    int file_timing_exists = validate_if_file_exists(filename_timing);
    if ( file_timing_exists == 0 ) 
    {   
        fp_timing= fopen(filename_timing, "w" );
        init_timing_file_contents(fp_timing); 
    } else 
    {
        fp_timing = fopen(filename_timing, "a " );
    } 
    fprintf(stdout, "\n\n# RUNNING : \t%s %.2s %d %d \n", argv[0], init_type, nx, nt);

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

//  CALCULATION :: |C| and infinity norm using manual
    fprintf(stdout,"# RESULTS : complex manual calculation ... \n");
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C| for manual calculation ... \n", nx, ny);
    init_matrix_zero(C, nx, ny);
    print_matrix(C, nx, ny, fp_matrix);
    gettimeofday(&tv1, &tz);
#pragma omp parallel shared (A, B, C) private (omp_nt, ni, nj, nk, manual_row_norm, manual_inf_norm)
    {        
        omp_nt = omp_get_thread_num();
        int slice_start = (((omp_nt+1) - 1) * number_of_row_slices);
        int slice_end = (((omp_nt+1) * number_of_row_slices) - 1);
    
       //  for (ni=slice_start; ni<slice_end; ni++)         
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
    } // end parallel
#pragma omp critical
    {
        for (ni=0; ni<nx; ni++)
        {
            manual_row_norm =0.0;
            for (nj=0; nj<nx; nj++)
            {
                manual_row_norm += C[(ni*nx) +nj];
            }
            manual_inf_norm = ( manual_inf_norm < manual_row_norm ) ? manual_row_norm : manual_inf_norm;
        }
    } // end critcal
    gettimeofday(&tv2, &tz);
    manual_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# RESULTS : \n\t# |C| : <%d> x <%d> matrix computed values : MANUAL ... \n", nx, ny);
    print_matrix(C, nx, ny, fp_matrix);
    fprintf(fp_matrix, "\t# |C| : matrix infinity norm is %g, \tcalculated in %f seconds ... \n", manual_inf_norm, manual_elapsed);

//  CALCULATION :: |C| and infinity norm using DGEMM
    fprintf(stdout,"# RESULTS : BLAS/ATLAS calculation -\n");
    fprintf(fp_matrix, "# \n# Initialize matrix <%d> x <%d> |C|, redone for CBLAS/ATLAS ... \n", nx, ny);
    init_matrix_zero(C, nx, ny);
    print_matrix(C, nx, ny, fp_matrix);    
    multipy_ABC_cblas(A, B, C, nx, ny, ALPHA, BETA);
    for (ni=0; ni<nx; ni++)
    {
        dgemm_row_norm =0.0;
        for (nj=0; nj<nx; nj++)
        {
            dgemm_row_norm += C[(ni*nx) +nj];
        }
        dgemm_inf_norm = ( dgemm_inf_norm < dgemm_row_norm ) ? dgemm_row_norm : dgemm_inf_norm;
    }
    gettimeofday(&tv2, &tz);
    dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# RESULTS : \n\t# |C| : <%d> x <%d> matrix computed values : CBLAS/ATLAS ... \n", nx, ny);
    print_matrix(C, nx, ny, fp_matrix);
    fprintf(fp_matrix, "\t# |C| : matrix infinity norm is %g, \tcalculated in %f seconds ... \n", manual_inf_norm, manual_elapsed);

//  OUTPUT :: results to stdout & .dat file : matrix size || infinity norm || Time/manual || infinity norm || Time / dgemm
    fprintf(stdout,"# \t\t|Matrix|  |Threads|  Time/manual  Inf Norm/manual  Time/dgemm  Inf Norm/dgemm\n");
    fprintf(stdout,"# Results: \t%d \t%d \t%fs \t%g \t%fs \t%g \n", nx, nt, manual_elapsed, manual_inf_norm, dgemm_elapsed, dgemm_inf_norm);
    fprintf(fp_timing, "%d \t%d \t%fs \t%g \t%fs \t%g \n", nx, nt, manual_elapsed, manual_inf_norm, dgemm_elapsed, dgemm_inf_norm);
   	
//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(A);
    deallocate_matrix_memory(B);
    deallocate_matrix_memory(C);    
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
