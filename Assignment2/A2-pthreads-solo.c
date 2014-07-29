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
#include <pthread.h>
#include <string.h>
#include <cblas.h>

typedef struct 
{
    double *sliceB;
    double *sliceC;
    int num_of_rows;
    pthread_mutex_t *mutex;
}   matrix_slice;

double *allA;
double *allB;
double *allC;
double pthread_norm;
int nx, ny, nt;

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
    fprintf(fp, "# Program :\tA2-pthreads-solo \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "# --------------------------------------------------------------------------------------------------  \n ");
}

void initialize_data_file_contents (FILE *fp) 
{
    fprintf(fp, "# --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA2-pthreads-solo \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
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

// pthreads calculation
void *pthread_matrix_slice_multiply(void *arg)
{
    matrix_slice *slice = arg;
    int ni, nj;
    const double ALPHA = 1.0;
    const double BETA = 0.0;
    
//  compute |C|    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nx, slice->num_of_rows, nx, ALPHA, allA, nx, slice->sliceB, nx, BETA, slice->sliceC, nx);
//  compute row norm of each slice
    double slice_norm = 0.0;
    for (ni = 0; ni < slice->num_of_rows; ni++) 
    {
        double row_sum=0.;
        for (nj = 0; nj < nx; nj++)
        {
            row_sum += *(slice->sliceC + nj * nx + ni);
        }
        if  (row_sum > slice_norm)
        {
            slice_norm = row_sum;
        }
    }
    pthread_mutex_lock(slice->mutex);
    if (slice_norm > pthread_norm)
    {
        pthread_norm = slice_norm;
    }
    pthread_mutex_unlock(slice->mutex);
    pthread_exit(NULL);
}

nt main (int argc, char *argv[])
{

//  define variables
    struct timeval tv1, tv2;
    struct timezone tz;
    int increment_or_random = 0;    // increment enabled by default
    int max_num_args = 5;
    char filename_matrix[60];
    char filename_timing[60];
    int MAXN = 1000;
    int MAXT = 100;

    int nt, num_per_slice;
    pthread_t *working_thread;
    matrix_slice *slice;
    pthread_mutex_t *mutex;
    int ni = 0;

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
    nx = atoi(argv[2]);                                     
    if (nx > MAXN)
    {
        fprintf(stderr, "\nWARNING:\tMatrix size entered <nx> [%d]  too large, now set to [%d]. \n", nx, MAXN); 
        nx = MAXN;
    }       
// number of threads
    nt = atoi(argv[3]);                                     
    if (nt > MAXT)
    {
        fprintf(stderr, "\nWARNING:\tNumber of Threads entered <nt> [%d]  too large. \n\t\tUsing defaults <nx> [%d] & <nt> [%d]. \n", nt, MAXN, MAXT); 
        nx = MAXN;
        nt = MAXT;
    }    
    if (nx % nt)
    {
        fprintf(stderr, "\nWARNING:\tNumber of Threads entered <nt> [%d]  must divide evenly into <nx> [%d]. \n\t\tUsing defaults <nx> [%d] & <nt> [%d]. \n", nt, nx, MAXN, MAXT); 
        nx = MAXN;
        nt = MAXT;
    }    
//  number of columns    
    ny = nx;         
//  matrix file name .txt
    strncpy(filename_matrix, argv[4], 59);
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
    strncpy(filename_timing, argv[5], 59);
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
    fprintf(fp_matrix, "# \n# RUNNING :\t%s %.2s %d %d\n", argv[0], init_type, nx, nt);
    fprintf(stdout, "\n# RUNNING :\t%s %.2s %d %d \n", argv[0], init_type, nx, nt);

//  CREATE & INITIALIZE :: matrices allA & allB & allC and output results to matrix file for reference
    fprintf(stdout, "# ALLOCATE :\tmatrices |allA|, |allB| ... \n") ; 
    fprintf(fp_matrix, "\n# ALLOCATE :\tmatrices |allA|, |allB| ... \n") ; 
    double *allA = allocate_memory_matrix(nx, ny);
    double *allB = allocate_memory_matrix(nx, ny);
    double *allC = allocate_memory_matrix(nx, ny);
    fprintf(stdout,"# INITIALIZE :\t|allA| & |allB| ... \n");
    fprintf(fp_matrix, "# INITIALIZE :\t|allA| & |allB| ... \n");
    if (increment_or_random == 1) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |allA|  using random mod 10 ... \n", nx, ny);
        initialize_matrix_random(allA, nx, ny);
        print_matrix(allA, nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |allB| using random mod 10 ... \n", nx, ny);
        initialize_matrix_random(allB, nx, ny);
        print_matrix(allB, nx, ny, fp_matrix);
    } else if (increment_or_random == 0) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |allA| using incremental <column> value + 1 ... \n", nx, ny);
        initialize_matrix_increment(allA, nx, ny);
        print_matrix(allA, nx, ny, fp_matrix);            
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |allB| using incremental <column> value + 1 ... \n", nx, ny);
        initialize_matrix_increment(allB, nx, ny);
        print_matrix(allB, nx, ny, fp_matrix);            
    }
    
//  MANUAL execution : calculate |allC| and infinity norm for resulting |allC|
    fprintf(stdout,"# INITIALIZE :\t<%d> x <%d> matrix |allC| for Straight-forward IJK manual computation ... \n", nx, ny);
    fprintf(fp_matrix, "# INITIALIZE :\t<%d> x <%d> matrix |allC| for Straight-forward IJK manual computation ... \n", nx, ny) ; 
    initialize_matrix_zero(allC, nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double manual_norm = multipy_abc_manual(allA, allB, allC, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double manual_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6 ; 
    fprintf(fp_matrix, "# RESULTS :\tmanual Straight-forward IJK calculation ... \n") ; 
    fprintf(fp_matrix, "# \t\tComputed Matrix [%d] x [%d] |allC| ... \n", nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ;     
    fprintf(fp_matrix, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", manual_elapsed, manual_norm) ; 
    fprintf(stdout, "# RESULTS :\tmanual Straight-forward IJK calculation ... \n") ; 
    fprintf(stdout, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%.1f] ... \n", manual_elapsed, manual_norm) ; 


//  PTHREADS execution : calculate |allC| and infinity norm for resulting |allC|    
    fprintf(stdout,"# INITIALIZE :\t|allC| for BLAS/ATLAS computation ... \n");
    fprintf(fp_matrix, "# INITIALIZE :\t<%d> x <%d> matrix |allC| for BLAS/ATLAS computation ... \n", nx, ny) ; 
    initialize_matrix_zero(allC, nx, ny) ; 
    print_matrix(allC, nx, ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz);
    working_thread = malloc(nt * sizeof(pthread_t));
    slice = malloc(nt * sizeof(matrix_slice));
    mutex = malloc(sizeof(pthread_mutex_t));
    num_per_slice = nx / nt;
    for (ni = 0; ni < nt; ni++)
    {
        slice[ni].sliceB = allB + ni * num_per_slice;
        slice[ni].sliceC = allC + ni * num_per_slice;
        slice[ni].mutex = mutex;
        slice[ni].num_of_rows = (ni == nt - 1) ? nx-ni * num_per_slice : num_per_slice;
        pthread_create(&working_thread[ni], NULL, pthread_matrix_slice_multiply, (void *)&slice[ni]);
    }
    for (ni = 0; ni < nt; ni++)
    {
        pthread_join(working_thread[ni], NULL);
    }
    gettimeofday(&tv2, &tz);
    double pthread_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
    printf("elapsed time: %f \n", pthread_elapsed);   
    printf("row sum norm is %f \n", pthread_norm);  

//  OUTPUT :: results to stdout & .dat file : |Matrix| ||inf norm/manual || Time/manual || infinity norm || Time / dgemm
    fprintf(stdout,"# SUMMARY :\t|Matrix|  \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \n");
    fprintf(stdout,"# \t\t%d \t\t%f \t%.1f \t\t%f \t%.1f \n", nx, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm);
    fprintf(fp_timing, "%d \t%f \t%.1f \t%f \t%.1f \n", nx, manual_elapsed, manual_norm, dgemm_elapsed, dgemm_norm);
   	
//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(allA);
    deallocate_matrix_memory(allB);
    deallocate_matrix_memory(allC);    
    deallocate_matrix_memory(working_thread);
    deallocate_matrix_memory(slice);
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
