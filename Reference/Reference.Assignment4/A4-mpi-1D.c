/* 
* *******************************************************************************

    Paula Dwan
    13208660 : COMP40730 : Assignment 4

    MPI program computing norm of product of two nxn matrices
    1/     Granularity of the product : 
           (ii)      One step algorithim. No intermediate resulting matrix
    2/     Partitioning Scheme :
           (a)      Left matrix is horizonitally partitioned
    3/     Matrix norm to be computed :
           (b)      Maximum absolute row sum norm (aka infinity-norm)
                     calculate max value of sum of each row (or product of same) 
                     ||A|| infintity = max SUM | Aij |
                     where   max     :==> 0 <= i <=n
                                  SUM   :==>  j =0 to n-1      
                                  
* ********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>
#include <string.h>
#include <mpi.h>

void usage () 
{
    fprintf(stdout,"\nUSAGE : \t<program name> [<-r>|<-i>] [N] <matrix contents file>.txt <timing file>.dat \n");
    fprintf(stdout,"\nTO : \t\tCalculate |C| = |A| x |B| using MPI and also calculate infinity norm of |C|. \n");    
    fprintf(stdout,"\nWHERE :");
    fprintf(stdout,"\t1. \t<-r>\tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n");
    fprintf(stdout,"\t \t<-i>\tinitialize |A| & |B| _incrementally_ with <column> value and |C| with '0' \n");
    fprintf(stdout,"\t2. \t[N] \tmax size of each matrix, if invalid defaults to 1,000 \n");
    fprintf(stdout,"\t3. \t[P] \tnumber of processors, (i) less than [N] and (ii) [N] mod [P] = 0\n");
    fprintf(stdout,"\t4. \t<matrix contents file>.txt\n\t\tname of .txt file to store values of matrices |A| |B| & |C| \n");
    fprintf(stdout,"\t5. \t<timing .dat file> .dat \n\t\tname of .dat file to contain time to complete for each iteration \n \n");
    exit(0);
}

int validate_if_file_exists(char * fn)
{
    FILE *fp = fopen(fn, "r") ; 
    if ( fp !=  NULL )
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
        fprintf(stderr,"ERROR : \texiting - memory allocation failed for matrix.\n");
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
    for (ni=0; ni<(rows*cols); ni++)
    {
        fprintf(fp, "%g\t", l_matrix[ni]); 
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

void init_matrix_random(double *l_matrix, int rows, int cols)
{
    int ni, mod_rand=10;
    for (ni=0; ni<(rows*cols); ni++)
    {
        l_matrix[ni]= (double) ((rand()%mod_rand) + 1);
    }
}

void init_matrix_increment(double *l_matrix, int rows, int cols)
{
    int ni;    
    for (ni=0; ni<(rows*cols); ni++)
    {
        l_matrix[ni]= ((ni % rows) + 1);
    }
}

void init_matrix_zero(double *l_matrix, int rows, int cols)
{
    int ni;
    for (ni=0; ni<(rows*cols); ni++)
    {
        l_matrix[ni]=  0.0;
    }
}

void init_matrix_file (FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA4-mpi-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n ");
}

void init_data_file (FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA4-mpi-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n");
    fprintf(fp, "# |Matrix| \t|Threads| \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \tTime/mpi \tInf Norm/mpi\n# \n");
}

//  manual calculation
double multipy_abc_manual(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols)
{
    int ni, nj, nk ; 
    
    for (ni=0 ; ni<rows ; ni++)
    {
        for (nj=0 ; nj<cols ; nj++)
        {
            double sum = 0.0 ; 
            for (nk=0 ; nk<rows ; nk++)
            {
                sum+= (matrix_a[(ni*rows)+nk]) * (matrix_b[(nk*rows)+nj]) ; 
            }
            matrix_c[(ni*rows)+nj] = sum ; 
        }
    }
double row_norm =0.0, inf_norm =0.0 ; 
    for (ni=0 ; ni<rows ; ni++)
    {
        row_norm =0.0 ; 
        for (nj=0 ; nj<rows ; nj++)
        {
            row_norm += matrix_c[(ni*rows) +nj] ; 
        }
        inf_norm = ( inf_norm < row_norm ) ? row_norm : inf_norm ; 
    }
    return inf_norm ; 
}

//  dgemm calculation 
double multipy_abc_dgemm(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols)
{
    int ni, nj ; 
// m, n, k :    local integers indicating the size of the matrices for rows x columns :: A :  m x  k, B :  k x n, C:  m x n
//                 Here, m = n = k = rows = columns = <nx> = <ny> as supplied
    int lm = rows, ln = rows ; 
// la_offset, lb_offset, lc_offset :
//                 Leading dimension of matrix A, B or C respectively, or the number of elements between 
//                 successive rows for row-major storage or columns for column-major storage. 
    int la_offset = rows, lb_offset = cols, lc_offset = rows ; 
    int ALPHA=1.0 ; 
    int BETA=0.0 ; 
    
    cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, lm, ln, ln, ALPHA, matrix_a, la_offset, matrix_b, lb_offset, BETA, matrix_c, lc_offset) ;   

    double row_norm =0.0, inf_norm =0.0 ; 
    for (ni=0 ; ni<rows ; ni++)
    {
        row_norm =0.0 ; 
        for (nj=0 ; nj<rows ; nj++)
        {
            row_norm += matrix_c[(ni*rows) +nj] ; 
        }
        inf_norm = ( inf_norm < row_norm ) ? row_norm : inf_norm ; 
    }
    return inf_norm ; 
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main ( int argc, char *argv[] )
{

//  define variables
    struct timeval tv1, tv2;
    struct timezone tz;
    int increment_or_random = 0;    // random enabled by default
    int max_num_args = 6;
    char filename_matrix[50];
    char filename_timing[50];
    int MAXN = 1000;
    int MAXP = 10;
    int size, rank, ni, nj, nk;

//  MPI Initialization
    MPI_init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
          
//  CLI PARAMETERS :: validate and initialize
    if ( argc != max_num_args ) 
    {
        fprintf(stderr, "\nERROR: \t<number of arguments> %d : is invalid, less than <default> %d\n", argc, max_num_args);      
        MPI_Finalize();  
        usage();
    }  
//  random or increment initialization of matrices |A| and |B|
    char cell_value_set[3];
    strncpy(cell_value_set, argv[1], 2);
    cell_value_set[3] = '\0';
    if (strcmp(cell_value_set, "-i") == 0)
    {
        increment_or_random = 1 ;        
    } else if (strcmp(cell_value_set, "-r") == 0)
    {
        increment_or_random = 0 ; 
    } else 
    {
        fprintf(stderr, "\nERROR : \tNeither switch '-i' nor '-r' entered, exiting. \n") ;
        MPI_Finalize();  
        usage(); 
    }
//  matrix size
    int nx = atoi(argv[2]);                                     
    if ( nx > MAXN )
    {
        fprintf(stderr, "\nWARNING: \tMatrix size entered <nx> %d  too large, now set to %d \n", nx, MAXN); 
        nx = MAXN;
    }    
    int ny = nx;        
//  number of processors 
    int np = atoi(argv[3]);  
    if ( (nx % np) != 0)
    {
        fprintf(stderr, "\nWARNING : \t<np> [%d] : number of processors must divide evenly into <nx> [%d] : matrix size. \n\t\tUsing defaults <np> [%d] & <nx> [%d] \n", np, nx, MAXP, MAXN) ; 
        nx=MAXN ; 
        np=MAXP ;        
    }
    if (np > nx)
    {
        fprintf(stderr, "\nWARNING: \t<np> [%d] : number of processors must be less than <nx> [%d] : matrix size. \n\t\tUsing defaults <np> [%d] & <nx> [%d] \n", np, nx, MAXP, MAXN) ; 
        nx=MAXN ; 
        np=MAXP ; 
    }
//  matrix file name .txt
    strncpy(filename_matrix, argv[4], 49);
    filename_matrix[50] = '\0';
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
    int file_timing_exists = validate_if_file_exists(filename_timing);    
    if ( file_timing_exists == 0 ) 
    {   
        fp_timing= fopen(filename_timing, "wa" );
        init_timing_file_contents(fp_timing); 
    } else 
    {
        fp_timing = fopen(filename_timing, "a" );
    }
    fprintf(fp_matrix, "# \n# RUNNING : \t%s %s %d %d %s %s \n", argv[0], cell_value_set, nx, np);
    fprintf(stdout, "\n# RUNNING : \t%s %s %d %d %s %s \n", argv[0], cell_value_set, nx, np);

//  CREATE & INITIALIZE :: matrices A & B & C and output results to matrix file for reference
    if ( rank == 0 ) 
    {    
        fprintf(stdout, "CREATE MATRICES |A|, |B| and |C| ... \n") ; 
        fprintf(fp_matrix, "# CREATE MATRICES |A|, |B| and |C| ... \n") ; 
        double *A =  allocate_memory_matrix(nx, ny);
        double *B =  allocate_memory_matrix(nx, ny);
        double *C =  allocate_memory_matrix(nx, ny);
        fprintf(stdout,"# INITIALIZE MATRICES |A|, |B| and |C| ... \n");
        fprintf(fp_matrix, "# INITIALIZE MATRICES |A|, |B| and |C| ... \n") ; 
        if (increment_or_random == 0) 
        {
            fprintf(fp_matrix, "\n# Initialize results <%d> x <%d> |A| using random mod 10 ... \n", nx, ny);
            init_matrix_random(A, nx, ny);
            fprintf(fp_matrix, "\n# Initialize results <%d> x <%d> |B| using random mod 10 ... \n", nx, ny);
            init_matrix_random(B, nx, ny);
      } else if (increment_or_random == 1) 
        {
            fprintf(fp_matrix, "\# Initialize results <%d> x <%d> |A| using incremental <column> value + 1 ... \n", nx, ny);
            init_matrix_increment(A, nx, ny);
            fprintf(fp_matrix, "\# Initialize results <%d> x <%d> |B| using incremental <column> value + 1 ... \n", nx, ny);
            init_matrix_increment(B, nx, ny);
        }
        print_matrix(A, nx, ny, fp_matrix);
        print_matrix(B, nx, ny, fp_matrix);
    }

//  MANUAL execution : calculate |C| and infinity norm for resulting |C|
    fprintf(fp_matrix, "# Initialize matrix [%d] x [%d] |C| for Straight-forward IJK manual computation ... \n", nx, ny) ; 
    init_matrix_zero(C, nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double manual_norm = multipy_abc_manual(A, B, C, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double manual_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6 ; 
    fprintf(fp_matrix, "\n# RESULTS : manual Straight-forward IJK calculation ... \n") ; 
    fprintf(fp_matrix, "# Computed Matrix [%d] x [%d] |C| ... \n", nx, ny) ; 
    fprintf(fp_matrix, "# |C| : matrix calculated in %f seconds and has infinity norm of [%g]  ... \n", manual_elapsed, manual_norm) ; 
    fprintf(stdout, "# RESULTS :  manual Straight-forward IJK calculation ... \n") ; 
    fprintf(stdout, "# |C| : matrix calculated in %f seconds and has infinity norm of [%g]  ... \n", manual_elapsed, manual_norm) ; 

//  CALCULATION :: |C| using DGEMM
    fprintf(fp_matrix, "# Initialize matrix [%d] x [%d] |C| for BLAS/ATLAS computation ... \n", nx, ny) ; 
    init_matrix_zero(C, nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double dgemm_norm = multipy_abc_dgemm(A, B, C, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6 ; 
    fprintf(stdout, "# RESULTS : BLAS/ATLAS computation ...\n") ; 
    fprintf(stdout, "# Matrix |C| calculated in %f seconds and has infinity norm of [%g]  ... \n", dgemm_elapsed, dgemm_norm) ;     
    fprintf(fp_matrix, "\n# RESULTS : BLAS/ATLAS computation ...\n") ; 
    fprintf(fp_matrix, "\n# Computed Matrix [%d] x [%d] |C| ... \n", nx, ny) ; 
    fprintf(fp_matrix, "# Matrix |C| calculated in %f seconds and has infinity norm of [%g]  ... \n", dgemm_elapsed, dgemm_norm) ; 

//  CALCULATION :: |C| using MPI
    fprintf(fp_matrix, "# Initialize matrix [%d] x [%d] |C| for MPI calculation ... \n", nx, ny) ; 
    init_matrix_zero(C, nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ; 
    int nx_proc_segment = ( nx * ny ) / np;
    int nx_start = ( rank * size ) / np;
    int nx_end = ( ( rank + 1 ) * size ) / np;    
    gettimeofday(&tv1, &tz) ; 
    MPI_Bcast (B, nx*ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter (A, nx_proc_segment, MPI_DOUBLE, A+nx_start, nx_proc_segment, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (ni=nx_start; ni<nx_end; ni++)
    {
        for (nj=0; nj<ny; nj++)
        {
            for (nk=0; nk<nx; nk++)
            {
                C[(ni*nx)+nj] += (A[(ni*nx)+nk]) * (B[(nk*nx)+nj]);
            }
        }
    }
    MPI_Gather (C+nx_end, nx_proc_segment, MPI_DOUBLE, C, nx_proc_segment, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double mpi_norm =0.0,mpi_row_norm =0.0 ; 
    if (rank ==0 )
    {
        for (ni=istart ; ni<iend ; ni++)
        {
            mpi_row_norm += matrix_c[ni] ; 
        }
        mpi_norm = ( mpi_norm <mpi_row_norm ) ? mpi_row_norm  : mpi_norm; 
    }
    gettimeofday(&tv2, &tz) ; 
    double mpi_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6 ; 
    print_matrix(C, nx, ny, stdout) ;     
    fprintf(stdout, "# RESULTS : MPI computation ...\n") ; 
    fprintf(stdout, "# Matrix |C| calculated in %f seconds and has infinity norm of [%g]  ... \n", mpi_elapsed, mpi_norm) ;     
    fprintf(fp_matrix, "\n# RESULTS : Straight-forward IJK BLAS/ATLAS computation ...\n") ; 
    fprintf(fp_matrix, "\n# Computed Matrix [%d] x [%d] |C| ... \n", nx, ny) ; 
    print_matrix(C, nx, ny, fp_matrix) ; 
    fprintf(fp_matrix, "# Matrix |C| calculated in %f seconds and has infinity norm of [%g]  ... \n", mpi_elapsed, mpi_norm) ; 

//  OUTPUT :: results to stdout & .dat file : |Matrix| || |Processors| ||inf norm/manual || Time/manual || infinity norm || Time / dgemm
    fprintf(stdout,"\t|Matrix| \t|Threads| \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \tTime/mpi \tInf Norm/mpi\n");
    fprintf(stdout,"Results: \t%d \t \t%lfs \t%g \t%lfs \t%g \n", nx, np, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm,mpi_elapsed, mpi_norm);
    fprintf(fp_timing, "%d \t%lfs \t%g \t%lfs \t%g \n", nx, np, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm);
   	
//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(A);
    deallocate_matrix_memory(B);
    deallocate_matrix_memory(C);    
    fclose(fp_matrix);
    fclose(fp_timing);
    MPI_Finalize();  
    return 0;
}
