/* 
* *******************************************************************************

    Paula Dwan
    13208660 : COMP40730 : Assignment 4

    MPI program computing norm of product of two nxn matrices
    1/     Granularity of the product : 
           (ii) One step algorithim. No intermediate resulting matrix
    2/     Partitioning Scheme :
           (a)  Left matrix is horizonitally partitioned
    3/     Matrix norm to be computed :
           (b)  Maximum absolute row sum norm (aka infinity-norm)
                calculate max value of sum of each row (or product of same) 
                |A| infintity = max SUM | Aij |
                where   max :==>    0 <= i <= n
                        SUM :==>    0 <= j <= n-1      
                        
    This program A4-mpi-manual.c computes the time taked for manual evaulation of |allC| and 
    the infinity norm of |allC| using :
    (i)     Straight-forward IJK 
    (ii)    DGEMM computations.                        
    
    The program A4-mpi-solo.c calculates |allC| and the infinity norm of |allC| using :
    (i)     Straight-forward IJK 
    (ii)    MPI 
                                  
* ********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cblas.h>
#include <sys/time.h>

int nx, ny;

void usage () 
{
    fprintf(stdout,"\nUSAGE :\t<program name> [<-r>|<-i>] [N] <matrix contents file>.txt <timing file>.dat \n");
    fprintf(stdout,"\nTO :\tCalculate |allC| = |A| x |B| using MPI and also calculate infinity norm of |allC|. \n");    
    fprintf(stdout,"\nWHERE :\t1. \t<-r> \tinitialize |A| & |B| with _random_ numbers and |allC| with '0'. \n");
    fprintf(stdout," \t \t<-i> \tinitialize |A| & |B| _incrementally_ with <column> value and |allC| with '0'. \n");
    fprintf(stdout," \t2. \t[N] \tmax size of each matrix, if invalid defaults to 100. \n");
    fprintf(stdout," \t3. \t<matrix contents file>.txt\n \t \tname of .txt file to store values of matrices |A|, |B| & |allC| \n");
    fprintf(stdout," \t4. \t<timing .dat file> .dat \n \t \tname of .dat file to contain time to complete for each iteration \n\n");
    fprintf(stdout,"\tMANUAL\tStraight-forward IJK & DGEMM computations only.\n");
    fprintf(stdout,"\tSOLO\tStraight-forward IJK & MPI computations only. \n\n");
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

int validate_if_file_exists(char * fn)
{
    FILE *fp = fopen(fn, "r");
    if (fp !=  NULL)
    {
        fclose(fp);
        return 1;      // file exists
    }
    return 0;          // files does not exist
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
        l_matrix[ni]= (double) ((rand()%mod_rand) + 1);
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
    fprintf(fp, "# Program : \tA4-mpi-manual \n# where : \t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "# -------------------------------------------------------------------------------------------------- \n");
}

void initialize_data_file_contents (FILE *fp) 
{
    fprintf(fp, "# --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program : \tA4-mpi-manual \n# where : \t.dat contains timing data & .txt contains matrix values\n# \n");
    fprintf(fp, "# --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# |Matrix| \tTime/dgemm \tInf Norm/dgemm \tTime/manual \tInf Norm/manual\n# \n");
}

double multipy_abc_manual(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols)
{
    int ni, nj, nk;
    double row_norm =0.0, inf_norm =0.0;
    
    for (ni=0;ni<rows;ni++)
    {
        for (nj=0;nj<cols;nj++)
        {
            double sum = 0.0;
            for (nk=0;nk<rows;nk++)
            {
                sum+= (matrix_a[(ni*rows)+nk]) * (matrix_b[(nk*rows)+nj]);
            }
            matrix_c[(ni*rows)+nj] = sum;
        }
    }
    for (ni=0;ni<rows;ni++)
    {
        row_norm =0.0;
        for (nj=0;nj<rows;nj++)
        {
            row_norm += matrix_c[(ni*rows) +nj];
        }
        inf_norm = (inf_norm < row_norm) ? row_norm : inf_norm;
    }
    return inf_norm;
}

double multipy_abc_dgemm(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols)
{
    int ni, nj;
    double ALPHA=1.0, BETA=0.0; 
    double row_norm =0.0, inf_norm =0.0;
//  m, n, k : local integers indicating the size of the matrices for rows x columns :: A :  m x  k, B :  k x n, C:  m x n
//  Here, m = n = k = rows = columns = <nx> = <ny> as supplied
    int lm = rows, ln = rows;
//  la_offset, lb_offset, lc_offset : Leading dimension of matrix A, B or C respectively, or the number of elements between 
//  successive rows for row-major storage or columns for column-major storage. 
    int la_offset = rows, lb_offset = cols, lc_offset = rows;
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lm, ln, ln, ALPHA, matrix_a, la_offset, matrix_b, lb_offset, BETA, matrix_c, lc_offset);  
    for (ni=0;ni<rows;ni++)
    {
        row_norm =0.0;
        for (nj=0;nj<rows;nj++)
        {
            row_norm += matrix_c[(ni*rows) +nj];
        }
        inf_norm = (inf_norm < row_norm) ? row_norm : inf_norm;
    }
    return inf_norm;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main (int argc, char *argv[])
{

//  define variables
    const int max_num_args = 5;
    const int MAXN = 1000;
    int increment_or_random = 0;    // increment enabled by default
    char filename_matrix[50];
    char filename_timing[50];
    double *allA, *allB, *allC;
    struct timeval tv1, tv2;
    struct timezone tz;
         
//  CLI PARAMETERS :: validate and initialize
    if (argc != max_num_args) 
    {
        fprintf(stderr, "\nERROR: \t<number of arguments> [%d] : is invalid, less than <default> [%d].\n", argc, max_num_args);      
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
        fprintf(stderr, "\nERROR :\t'invalid entry : %s for '-i' or '-r'", init_type);        
        usage();
    }

//  matrix size
    nx = atoi(argv[2]);                                     
    if (nx > MAXN)
    {
        fprintf(stderr, "\nWARNING :\tMatrix size entered <nx> [%d]  too large, now set to [%d]. \n", nx, MAXN); 
        nx = MAXN;
    }    
    ny = nx;        
//  matrix file name .txt
    strncpy(filename_matrix, argv[3], 49);
    filename_matrix[50] = '\0';
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
    strncpy(filename_timing, argv[4], 49);
    filename_timing[50] = '\0';    
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

//  ALLOCATE & INITIALIZE :: matrices and output results to matrix file for reference
    fprintf(stdout,"# ALLOCATE :\t|allA| and |allB| ... \n");
    fprintf(fp_matrix, "# ALLOCATE :\t|allA| and |allB| ... \n");
    allA = allocate_memory_matrix(nx, ny);           // needed for manual and dgemm
    allB = allocate_memory_matrix(nx, ny);
    fprintf(stdout,"# INITIALIZE :\tmatrices ... \n");
    fprintf(fp_matrix,"# INITIALIZE :\tmatrices ... \n");
    if (increment_or_random == 1) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |allA| using random mod 10 ... \n", nx, ny);
        initialize_matrix_random(allA, nx, ny);
        print_matrix(allA, nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t<%d> x <%d> matrix |allB| using random mod 10 ... \n", nx, ny);
        initialize_matrix_random(allB, nx, ny);
        print_matrix(allB, nx, ny, fp_matrix);
    } else if (increment_or_random == 0) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |allA| using incremental <column> value + 1... \n", nx, ny);
        initialize_matrix_increment(allA, nx, ny);
        print_matrix(allA, nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |allB| using incremental <column> value + 1 ... \n", nx, ny);
        initialize_matrix_increment(allB, nx, ny);
        print_matrix(allB, nx, ny, fp_matrix);
    }

//  MANUAL execution : calculate |allC| and infinity norm for resulting |allC|
        fprintf(stdout,"# INITIALIZE :\t<%d> x <%d> matrix |allC| for Straight-forward IJK manual computation ... \n", nx, ny);
        fprintf(fp_matrix, "# INITIALIZE :\t<%d> x <%d> matrix |allC| for Straight-forward IJK manual computation ... \n", nx, ny);
        initialize_matrix_zero(allC, nx, ny);
        print_matrix(allC, nx, ny, fp_matrix);
        gettimeofday(&tv1, &tz);
        double manual_norm = multipy_abc_manual(allA, allB, allC, nx, ny);
        gettimeofday(&tv2, &tz);
        double manual_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
        fprintf(stdout, "# RESULTS : \tStraight-forward IJK computation ...\n");
        fprintf(stdout, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", manual_elapsed, manual_norm);    
        fprintf(fp_matrix, "# RESULTS : \tStraight-forward IJK computation ...\n");
        fprintf(fp_matrix, "# \t\tComputed Matrix [%d] x [%d] |allC| ... \n", nx, ny);
        print_matrix(allC, nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n\n", manual_elapsed, manual_norm);

//  DGEMM execution : calculate |allC| and infinity norm for resulting |allC|
        fprintf(stdout,"# INITIALIZE :\t<%d> x <%d> matrix |allC| for BLAS/ATLAS computation ... \n", nx, ny);
        fprintf(fp_matrix, "# INITIALIZE :\t<%d> x <%d> matrix |allC| for BLAS/ATLAS computation ... \n", nx, ny);
        initialize_matrix_zero(allC, nx, ny);
        print_matrix(allC, nx, ny, fp_matrix);
        gettimeofday(&tv1, &tz);
        double dgemm_norm = multipy_abc_dgemm(allA, allB, allC, nx, ny);
        gettimeofday(&tv2, &tz);
        double dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
        fprintf(stdout, "# RESULTS : \tDGEMM IJK computation ...\n");
        fprintf(stdout, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", dgemm_elapsed, dgemm_norm);    
        fprintf(fp_matrix, "# RESULTS : \tDGEMM IJK computation ...\n");
        fprintf(fp_matrix, "# \t\tComputed Matrix [%d] x [%d] |allC| ... \n", nx, ny);
        print_matrix(allC, nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n\n", dgemm_elapsed, dgemm_norm);

//  OUTPUT :: results to stdout & .dat file
        fprintf(stdout,"# SUMMARY:\t|Matrix|  Time/dgemm Inf Norm/dgemm Time/manual Inf Norm/manual\n");
        fprintf(stdout,"\t\t%d \t%f \t%.1f \t%f \t%.1f f\n", nx, dgemm_elapsed, dgemm_norm, manual_elapsed, manual_norm);
        fprintf(fp_timing, "%d \t%f \t%.1f \t%f \t%.1f  \n", nx, dgemm_elapsed, dgemm_norm, manual_elapsed, manual_norm);

//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(allA); 
    deallocate_matrix_memory(allB);
    deallocate_matrix_memory(allC); 
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
