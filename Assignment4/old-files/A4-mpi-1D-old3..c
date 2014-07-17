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
                     ||allA|| infintity = max SUM | Aij |
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
    fprintf(stdout,"\nTO : \t \tCalculate |allC| = |allA| x |allB| using MPI and also calculate infinity norm of |allC|. \n");    
    fprintf(stdout,"\nWHERE :");
    fprintf(stdout," \t1. \t<-r> \tinitialize |allA| & |allB| with _random_ numbers and |allC| with '0' \n");
    fprintf(stdout," \t \t<-i> \tinitialize |allA| & |allB| _incrementally_ with <column> value and |allC| with '0' \n");
    fprintf(stdout," \t2. \t[N] \tmax size of each matrix, if invalid defaults to 1,000 \n");
    fprintf(stdout," \t3. \t<matrix contents file>.txt\n \t \tname of .txt file to store values of matrices |allA| |allB| & |allC| \n");
    fprintf(stdout," \t4. \t<timing .dat file> .dat \n \t \tname of .dat file to contain time to complete for each iteration \n \n");
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
        fprintf(fp, "%g \t", l_matrix[ni]); 
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

void init_matrix_file_contents (FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program : \tA4-mpi-1D \n# where : \t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n ");
}

void init_data_file_contents (FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program : \tA4-mpi-1D \n# where : \t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
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
        inf_norm = (inf_norm < row_norm) ? row_norm : inf_norm ; 
    }
    return inf_norm ; 
}

//  dgemm calculation 
double multipy_abc_dgemm(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols)
{
    int ni, nj ; 
// m, n, k :    local integers indicating the size of the matrices for rows x columns :: allA :  m x  k, allB :  k x n, allC:  m x n
//                 Here, m = n = k = rows = columns = <nx> = <ny> as supplied
    int lm = rows, ln = rows ; 
// la_offset, lb_offset, lc_offset :
//                 Leading dimension of matrix allA, allB or allC respectively, or the number of elements between 
//                 successive rows for row-major storage or columns for column-major storage. 
    int la_offset = rows, lb_offset = cols, lc_offset = rows ; 
    int ALPHA=1.0 ; 
    int BETA=0.0 ; 
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lm, ln, ln, ALPHA, matrix_a, la_offset, matrix_b, lb_offset, BETA, matrix_c, lc_offset) ;   

    double row_norm =0.0, inf_norm =0.0 ; 
    for (ni=0 ; ni<rows ; ni++)
    {
        row_norm =0.0 ; 
        for (nj=0 ; nj<rows ; nj++)
        {
            row_norm += matrix_c[(ni*rows) +nj] ; 
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
    int increment_or_random = 0;    // random enabled by default
    int max_num_args = 5;
    char filename_matrix[50];
    char filename_timing[50];
    int MAXN = 100;
    int MAXP = 10;
    int np, rank, ni, nj, nk;
    double sum;

//  MPI Initialization
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
          
//  CLI PARAMETERS :: validate and initialize
    if (argc != max_num_args) 
    {
        fprintf(stderr, "\nERROR: \t<number of arguments> %d : is invalid, less than <default> %d\n", argc, max_num_args);      
        MPI_Finalize();  
        usage();
    }  
//  random or increment initialization of matrices |allA| and |allB|
    char init_type[3];
    strncpy(init_type, argv[1], 2);
    init_type[3] = '\0';
    if (strcmp(init_type, "-i") == 0)
    {
        increment_or_random = 1 ;        
    }
    if (strcmp(init_type, "-r") == 0)
    {
        increment_or_random = 0 ; 
    } else  if (strcmp(init_type, "-r") == 1)
    {
        fprintf(stderr, "\nERROR : \tNeither switch '-i' nor '-r' entered. Exiting. \n") ;
        MPI_Finalize();  
        usage(); 
    }
//  matrix size
    int nx = atoi(argv[2]);                                     
    if (nx > MAXN)
    {
        fprintf(stderr, "\nWARNING: \tMatrix size entered <nx> [%d]  too large, now set to [%d]. \n", nx, MAXN); 
        nx = MAXN;
    }    
    int ny = nx;        
//  validate number of processors 
    if ((nx % np) != 0)
    {
        fprintf(stderr, "\nWARNING : \t<np> [%d] : number of processors must divide evenly into <nx> [%d] : matrix size. Existing. \n") ; 
        MPI_Finalize();  
        usage(); 
    }
    if (np > nx)
    {
        fprintf(stderr, "\nWARNING: \t<np> [%d] : number of processors must be less than <nx> [%d] : matrix size. Exiting. \n") ; 
        MPI_Finalize();  
        usage(); 
    }
//  matrix file name .txt
    strncpy(filename_matrix, argv[3], 49);
    filename_matrix[50] = '\0';
    FILE *fp_matrix;
    int file_matrix_exists = validate_if_file_exists(filename_matrix);
    if (file_matrix_exists == 0) 
    {   
        fp_matrix= fopen(filename_matrix, "wa");
        init_matrix_file_contents(fp_matrix); 
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
        init_data_file_contents(fp_timing); 
    } else 
    {
        fp_timing = fopen(filename_timing, "a");
    }
    fprintf(fp_matrix, "# \n# RUNNING : \t%s %.2s %d %d \n", argv[0], init_type, nx, np);
    fprintf(stdout, "\n# RUNNING : \t%s %.2s %d %d \n", argv[0], init_type, nx, np);

//  CREATE & INITIALIZE :: matrices allA & allB & allC and output results to matrix file for reference
    fprintf(stdout, "# CREATE MATRICES |allA|, |allB| and |allC| ... \n") ; 
    fprintf(fp_matrix, "\n\n# CREATE MATRICES |allA|, |allB| and |allC| ... \n") ; 
    double *allA = allocate_memory_matrix(nx, ny);
    double *allB = allocate_memory_matrix(nx, ny);
    double *allC = allocate_memory_matrix(nx, ny);
    fprintf(stdout,"# INITIALIZE : \t|allA| & |allB| ... \n");
    if (increment_or_random == 0) 
    {
        fprintf(fp_matrix, "# INITIALIZE : \t<%d> x <%d> matrix |allA|  using random mod 10 ... \n", nx, ny);
        init_matrix_random(allA, nx, ny);
        print_matrix(allA, nx, ny, fp_matrix);
        fprintf(fp_matrix, "# INITIALIZE : \t<%d> x <%d> matrix allB| using random mod 10 ... \n", nx, ny);
        init_matrix_random(allB, nx, ny);
        print_matrix(allB, nx, ny, fp_matrix);
    } else if (increment_or_random == 1) 
    {
        fprintf(fp_matrix, "# INITIALIZE : \t<%d> x <%d> matrix |allA| using incremental <column> value + 1 ... \n", nx, ny);
        init_matrix_increment(allA, nx, ny);
        print_matrix(allA, nx, ny, fp_matrix);            
        fprintf(fp_matrix, "# INITIALIZE : \t<%d> x <%d> matrix allB| using incremental <column> value + 1 ... \n", nx, ny);
        init_matrix_increment(allB, nx, ny);
        print_matrix(allB, nx, ny, fp_matrix);            
    }

//  MANUAL execution : calculate |allC| and infinity norm for resulting |allC|
    fprintf(stdout,"# INITIALIZE : \t|allC| ... \n");
    fprintf(fp_matrix, "# INITIALIZE : \t<%d> x <%d> matrix |allC| for Straight-forward IJK manual computation ... \n", nx, ny) ; 
    init_matrix_zero(allC, nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double manual_norm = multipy_abc_manual(allA, allB, allC, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double manual_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6 ; 
    fprintf(fp_matrix, "# RESULTS : \tmanual Straight-forward IJK calculation ... \n") ; 
    fprintf(fp_matrix, "# \t\tComputed Matrix [%d] x [%d] |allC| ... \n", nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ;     
    fprintf(fp_matrix, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n", manual_elapsed, manual_norm) ; 
    fprintf(stdout, "# RESULTS : \tmanual Straight-forward IJK calculation ... \n") ; 
    fprintf(stdout, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n", manual_elapsed, manual_norm) ; 

//  CALCULATION :: |allC| using DGEMM
    fprintf(stdout,"# INITIALIZE : \t|allC| ... \n");
    fprintf(fp_matrix, "# INITIALIZE : \t<%d> x <%d> matrix |allC| for BLAS/ATLAS computation ... \n", nx, ny) ; 
    init_matrix_zero(allC, nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double dgemm_norm = multipy_abc_dgemm(allA, allB, allC, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6 ; 
    fprintf(stdout, "# RESULTS : \tBLAS/ATLAS computation ...\n") ; 
    fprintf(stdout, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n", dgemm_elapsed, dgemm_norm) ;     
    fprintf(fp_matrix, "\n# RESULTS : \tBLAS/ATLAS computation ...\n") ; 
    fprintf(fp_matrix, "\n# \t\tComputed Matrix [%d] x [%d] |allC| ... \n", nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ; 
    fprintf(fp_matrix, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n", dgemm_elapsed, dgemm_norm) ; 

//  CALCULATION :: |allC| using MPI
    fprintf(stdout,"# INITIALIZE : \t|segmentA|, |segmentB|, |segmentC| and |allC| ... \n");
    fprintf(fp_matrix, "# INITIALIZE : \t|segmentA|, |segmentB|, |segmentC| and |allC| ... \n");
    int my_nx = nx / np;
    double *segmentA =  allocate_memory_matrix(my_nx, ny);
    double *segmentB =  allocate_memory_matrix(my_nx, ny);
    double *segmentC =  allocate_memory_matrix(my_nx, ny);
    fprintf(fp_matrix, "# \t<%d> x <%d> matrix  |segmentA| for MPI calculation ... \n", my_nx, ny) ; 
    for (ni=0; ni<my_nx; ni++) 
    {
        segmentA[ni] = allA[ni];
        fprintf(stdout, "DEBUG : \tsegmentA[ni] = allA[ni] : segmentA[%d] = allA[%d] : %g = %g \n", ni, ni, segmentA[ni], allA[ni]);
    }
    print_matrix(segmentA, my_nx, ny, fp_matrix) ; 
    fprintf(fp_matrix, "# \t<%d> x <%d> matrix  |segmentB| for MPI calculation ... \n", my_nx, ny) ; 
    for (ni=0; ni<my_nx; ni++) 
    {
        segmentB[ni] = allB[ni];
        fprintf(stdout, "DEBUG : \tsegmentB[ni] = allB[ni] : segmentB[%d] = allB[%d] : %g = %g \n", ni, ni, segmentB[ni], allB[ni]);
    } 
    print_matrix(segmentB, my_nx, ny, fp_matrix) ; 
    fprintf(fp_matrix, "# \t<%d> x <%d> matrix  |segmentC| for MPI calculation ... \n", my_nx, ny) ; 
    init_matrix_zero(segmentC, my_nx, ny) ; 
    print_matrix(segmentC, my_nx, ny, fp_matrix) ; 
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        gettimeofday(&tv1, &tz) ; 
    }
    for (ni=0; ni<np; ni++) {
        MPI_Gather(segmentB, nx*my_nx, MPI_DOUBLE, allB, nx*my_nx, MPI_DOUBLE, ni, MPI_COMM_WORLD);
    }
    for (ni=0; ni<my_nx; ni++) 
    {
        for (nj=0; nj<nx; nj++) 
        {
            sum= 0.0;
            for (nk<0; nk<nx; nk++) 
            {
                sum += segmentA[ni*nx + nk] * allB[nk*nx + nj];
            }
            segmentC [ni*nx + nj] = sum;
        }
    }
    deallocate_matrix_memory(allB);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)  
    {
    fprintf(stdout, "DEBUG : \tmatrix c now calculated \n");
    }
    if (rank == 0)  
    {
        fprintf(fp_matrix, "# INITIALIZE : \t<%d> x <%d> matrix  |allC| for MPI calculation ... \n", nx, ny) ; 
        init_matrix_zero(allC, nx, ny) ; 
        print_matrix(allC, nx, ny, fp_matrix) ; 
    }  
    MPI_Gather (segmentC, nx*my_nx, MPI_DOUBLE, allC, nx*my_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double mpi_norm =0.0, mpi_row_norm =0.0 , mpi_elapsed; 
    if (rank ==0)
    {
        for (ni=0 ; ni<nx ; ni++)
        {
            mpi_row_norm =0.0 ; 
            for (nj=0 ; nj<nx ; nj++)
            {
                mpi_row_norm += allC[(ni*nx) +nj] ; 
            }
            mpi_norm = (mpi_norm < mpi_row_norm) ? mpi_row_norm : mpi_norm ; 
        }
        gettimeofday(&tv2, &tz) ; 
        mpi_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6 ; 
    }
    fprintf(stdout, "# RESULTS : \tMPI computation ...\n") ; 
    fprintf(stdout, "# \t\t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n", mpi_elapsed, mpi_norm) ;     
    fprintf(fp_matrix, "\n# RESULTS : \tMPI computation ...\n") ; 
    fprintf(fp_matrix, "\n# \t\tComputed Matrix [%d] x [%d] |allC| ... \n", nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ; 
    fprintf(fp_matrix, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n", mpi_elapsed, mpi_norm) ; 

//  OUTPUT :: results to stdout & .dat file : |Matrix| || |Processors| ||inf norm/manual || Time/manual || infinity norm || Time / dgemm
    fprintf(stdout,"\t\t|Matrix|  |Processors| Time/manual Inf Norm/manual Time/dgemm \tInf Norm/dgemm Time/mpi Inf Norm/mpi\n");
    fprintf(stdout,"Results: \t%d \t%d \t%f \t%g \t%f \t%g \t%f \t%g \n", nx, np, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm, mpi_elapsed, mpi_norm);
    fprintf(fp_timing, "%d \t%d \t%f \t%g \t%f \t%g \t%f \t%g \n\n", nx, np, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm, mpi_elapsed, mpi_norm);
   	
//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    if (rank == 0) 
    {
        deallocate_matrix_memory(allC);    
    }
    MPI_Finalize();  
    // deallocate_matrix_memory(allA);
    deallocate_matrix_memory(segmentA);
    deallocate_matrix_memory(segmentB);
    deallocate_matrix_memory(segmentC);    
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
