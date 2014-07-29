/* 
* **********************************************************************************

    Paula Dwan
    13208660 : COMP40730 : Assignment 2

    Parallel pthreads program computing norm of product of two nxn matrices
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

* **********************************************************************************    
*/          

# include <sys/time.h>
# include <stdio.h>
# include <stdlib.h>
# include <pthread.h>
# include <cblas.h>
# include <math.h>
# include <string.h>

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
typedef struct 
{ 
     double *my_sliceA ; 
     double *my_sliceC ; 
     double *my_allB ; 
     double *my_norm_summation ; 
     double my_norm_sumC ; 
     pthread_mutex_t *mutex ; 
     int my_no_of_slices ; 
     int my_nx ; 
     int my_np ; 
} norm_sum_t ; 

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void usage (int maxN, int maxP) {
    fprintf(stdout, "\nUSAGE : \t<program name> [<-r>|<-i>] [N] [P] <matrix contents file>.txt <timing file>.dat \n") ; 
    fprintf(stdout, "\nTO : \t\tCalculate |C| = |A| x |B| using pthreads and also calculate infinity norm of |C| \n") ;    
    fprintf(stdout, "\nWHERE :") ; 
    fprintf(stdout, "\t1. \t<-r>\tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n") ; 
    fprintf(stdout, "\t   \t<-i>\tinitialize |A| & |B| _incrementally_ with <column> value and |C| with '0' \n") ; 
    fprintf(stdout, "\t2. \t[N] \tmatrix size, if invalid set to default Matrix Size [%d] \n", maxN) ; 
    fprintf(stdout, "\t3. \t[P] \tnumber of threads, i.e.: processors when (i) [P] less than [N] and (ii) [N] mod [P] = 0 \n\t\tIf invalid, set to defaults of Matrix Size [%d] and number of Threads [%d].  \n", maxN , maxP) ; 
    fprintf(stdout, "\t4. \t<matrix contents file>.txt\n\t\tname of .txt file to store values of matrices |A| |B| & |C| \n") ; 
    fprintf(stdout, "\t5. \t<timing .dat file> .dat \n\t\tname of timing data file to containing calculation time for each iteration \n\n") ; 
    exit(0) ; 
}

double* allocate_memory_matrix (int rows_cols)
{
    double* l_matrix ;  
    l_matrix = (double *) malloc(rows_cols*sizeof (double)) ;  
    if (!(l_matrix))
    {
        fprintf(stderr, "ERROR : \texiting - memory allocation failed for matrix.\n") ; 
        exit(1) ; 
    }
    return l_matrix ;  
} 

void deallocate_matrix_memory(double *matrix1d)
{
    free(matrix1d) ; 
}

void print_matrix(double *matrix1d, int rows_cols, FILE *fp)
{
    int ni ; 
    int count =1 ; 
    for (ni=0 ; ni<rows_cols ; ni++)
    {
        fprintf(fp, "%g\t", matrix1d[ni]) ; 
        if (count == sqrt(rows_cols) ) 
        {
            fprintf(fp, "\n") ; 
            count =1 ; 
        } else 
        {
            count ++ ;    
        }
    }
}

void init_matrix_random(double *matrix1d, int rows_cols)
{
    int ni, mod_rand=10 ; 
    for (ni=0 ; ni<(rows_cols) ; ni++)
    {
        matrix1d[ni]= (double) ((rand()%mod_rand) + 1) ; 
    }
}

void init_matrix_increment(double *matrix1d, int rows_cols)
{
    int ni ;    
    int n_rows_cols= sqrt(rows_cols) ; 
    for (ni=0 ; ni<(rows_cols) ; ni++)
    {
        matrix1d[ni]= ( ( ni % n_rows_cols ) + 1) ; 
    }
}

void init_matrix_zero(double *matrix1d, int rows_cols)
{
    int ni ; 
    for (ni=0 ; ni<(rows_cols) ; ni++)
    {
        matrix1d[ni]=  0.0 ; 
    }
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

void init_matrix_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n") ; 
    fprintf(fp, "# Program :\tA2-pthreads-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n") ; 
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n") ; 
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n") ; 
}

void init_timing_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n") ; 
    fprintf(fp, "# Program :\tA2-pthreads-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n") ; 
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n") ; 
    fprintf(fp, "# |Matrix| \t|Threads| \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \tTime/dgemm \tInf Norm/dgemm \n# \n") ; 
}

double multipy_abc_cblas(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols)
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

//  function to calculate infinity norm and provide values for |C| : manual
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

//  function to calculate infinity norm and provide values for |C|
void *norm_row_calculation(void *args)
{
     norm_sum_t *norm_sum_data ; 
     int i ;  
     int ALPHA=1.0 ; 
     int BETA=0.0 ; 
     
     norm_sum_data = args ; 
     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, norm_sum_data->my_nx, norm_sum_data->my_nx, norm_sum_data->my_nx, ALPHA, norm_sum_data->my_sliceA, norm_sum_data->my_nx, norm_sum_data->my_allB, norm_sum_data->my_nx, BETA, norm_sum_data->my_sliceC, norm_sum_data->my_nx) ; 
     fprintf(stdout, "# Evaluate sum of local |sliceC| using <mutex>...\n") ; 
     pthread_mutex_lock(norm_sum_data->mutex) ; 
     fprintf(stdout, "# DEBUG : \tnorm_sum_data->my_no_of_slices = [%d] \n ", (int)norm_sum_data->my_no_of_slices) ; 
     fprintf(stdout, "# DEBUG : \tvalue \tmax \n") ; 
     int iteration = (int) (norm_sum_data->my_no_of_slices / norm_sum_data->my_nx) ; 
     int counter = 0 ; 
     
     for(i=0 ; i< iteration ; i++) 
     {
          norm_sum_data->my_norm_summation[counter] += norm_sum_data->my_sliceC[i] ; 
          if ((counter % (int) norm_sum_data->my_nx) > 0)
          { 
               counter ++ ; 
          } else 
          {
               counter =0 ; 
          }
          if (norm_sum_data->my_norm_sumC < norm_sum_data->my_norm_summation[i]) 
          {
               norm_sum_data->my_norm_sumC = norm_sum_data->my_norm_summation[i] ; 
          }
     fprintf(stdout, "# DEBUG : \t%g \t %g \n",  (double)norm_sum_data->my_norm_summation[i], norm_sum_data->my_norm_sumC) ; 
     }  
     
     fprintf(stdout, "# DEBUG :  Max absolute value in local |sliceC| is : %g\n", norm_sum_data->my_norm_sumC) ;     
     pthread_mutex_unlock(norm_sum_data->mutex) ; 
     pthread_exit(NULL) ; 
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// start of main
int main ( int argc, char *argv[] )
{
//  define variables
    struct timeval tv1, tv2 ; 
    struct timezone tz ; 
    const int MAXPTHREADS = 100 ; 
    const int MAXN = 1000 ; 
    int i = 0 ; 
    int increment_or_random = 5 ; 
    int max_num_args = 6 ; 
    char filename_matrix[50] ; 
    char filename_timing[50] ; 
              
//  start of program

    if ( argc != max_num_args )
    {
        usage(MAXN, MAXPTHREADS) ; 
    } 
    
    char cell_value_type[3] ; 
    strncpy(cell_value_type, argv[1], 2) ; 
    cell_value_type[3] = '\0' ; 
    if (strcmp(cell_value_type, "-i") == 0)
    {
        increment_or_random = 1 ;        
    } else if (strcmp(cell_value_type, "-r") == 0)
    {
        increment_or_random = 0 ; 
    }
    int nx = atoi(argv[2]) ;                                     
    int ny = nx ;      
    int np = atoi(argv[3]) ;             
    if ( (nx % np) != 0)
    {
        fprintf(stderr, "\nWARNING : \t<np> [%d] : number of threads must divide evenly into <nx> [%d] : matrix size. \n\t\tUsing defaults <np> [%d] & <nt> [%d] \n", np, nx, MAXPTHREADS, MAXN) ; 
        nx=MAXN ; 
        np= MAXPTHREADS ;        
    }
    if (np > nx)
    {
        fprintf(stderr, "\nWARNING: \t<np> [%d] : number of threads must be less <nx> [%d] : matrix size. \n\t\tUsing defaults <np> [%d] & <nt> [%d] \n", np, nx, MAXPTHREADS, MAXN) ; 
        nx=MAXN ; 
        np= MAXPTHREADS ; 
    }
    
//  matrix file name .txt
    strncpy(filename_matrix, argv[4], 49) ; 
    filename_matrix[50] = '\0' ; 
    FILE *fp_matrix ; 
    int file_matrix_exists = validate_if_file_exists(filename_matrix) ; 
    if ( file_matrix_exists == 0 ) 
    {   
        fp_matrix= fopen(filename_matrix, "wa" ) ; 
        init_matrix_file_contents(fp_matrix) ; 
    } else 
    {
        fp_matrix = fopen(filename_matrix, "a+" ) ; 
    } 
//  data file name .dat
    strncpy(filename_timing, argv[5], 49) ; 
    filename_timing[50] = '\0' ;    
    FILE *fp_timing ; 
    int file_timing_exists = validate_if_file_exists(filename_timing) ; 
    if ( file_timing_exists == 0 ) 
    {   
        fp_timing= fopen(filename_timing, "wa" ) ; 
        init_timing_file_contents(fp_timing) ; 
    } else 
    {
        fp_timing = fopen(filename_timing, "a+" ) ; 
    } 
    fprintf(fp_matrix, "\n# RUNNING :\t%s %s %d %d \n", argv[0], cell_value_type, nx, np) ; 
    fprintf(stdout, "\n RUNNING : \t%s %s %d %d \n", argv[0], cell_value_type, nx, np) ; 
 
//  Calculate slice information using number of processors <np> and matrix size <nx>
    int no_of_slices = (nx*np) ; 
    int slice_size = (nx/np) ; 
    fprintf(fp_matrix, "# \tNo of slices = <matrix rows> x <number of processes> = [%d] x [%d] = [%d] ...\n", nx, np, no_of_slices) ; 
    fprintf(fp_matrix, "# \tSlice size       = <matrix rows> ÷ <number of processes> = [%d] ÷ [%d] = [%d] ...\n", nx, np, slice_size) ; 
    fprintf(stdout, "\tNo of slices = <matrix rows> x <number of processes> = [%d] x [%d] = [%d] ...\n", nx, np, no_of_slices) ; 
    fprintf(stdout, "\tSlice size       = <matrix rows> ÷ <number of processes> = [%d] ÷ [%d] = [%d] ...\n", nx, np, slice_size) ; 
 
// Allocate memory for matrices allA, allB allC, sliceA and sliceC
     fprintf(stdout, "CREATE MATRICES |allA|, |allB|, |allC|, |sliceA| & |sliceC| ... \n") ; 
     fprintf(fp_matrix, "# CREATE MATRICES |allA|, |allB|, |allC|, |sliceA| & |sliceC| ... \n") ; 
     double *allA = allocate_memory_matrix(nx*ny) ; 
     double *allB = allocate_memory_matrix(nx*ny) ; 
     double *allC = allocate_memory_matrix(nx*ny) ; 
     double *sliceA = allocate_memory_matrix(no_of_slices) ; 
     double *sliceC = allocate_memory_matrix(no_of_slices) ; 
 
 //  Initialize matrices A & B & C
    fprintf(stdout, "INITIALIZE MATRICES |allA|, |allB|, |allC|, |sliceA|  & |sliceC| ... \n") ; 
    fprintf(fp_matrix, "# INITIALIZE MATRICES |allA|, |allB|, |allC|, |sliceA|  & |sliceC| ... \n") ; 
    if (increment_or_random == 0) {
        fprintf(fp_matrix, "# Initialize : matrix [%d] x [%d] |allA| using <column> values ... \n", nx, ny) ; 
        init_matrix_random(allA, nx*ny) ; 
        fprintf(fp_matrix, "# Initialize : matrix [%d] x [%d] |allB| using <column> values ... \n", nx, ny) ; 
        init_matrix_random(allB, nx*ny) ; 
    } else if (increment_or_random == 1) 
    {
        fprintf(fp_matrix, "# Initialize : matrix [%d] x [%d] |allA| using random numbers 1 to 10 ... \n", nx, ny) ; 
        init_matrix_increment(allA, nx*ny) ; 
        fprintf(fp_matrix, "# Initialize : matrix [%d] x [%d] |allB|  using random numbers 1 to 10 ... \n", nx, ny) ; 
        init_matrix_increment(allB, nx*ny) ; 
    }
    print_matrix(allA, nx*ny, fp_matrix) ; 
    print_matrix(allB, nx*ny, fp_matrix) ; 
    fprintf(fp_matrix, "# Initialize : matrix [%d] x [%d] |sliceA| with Zero ... \n", 1, no_of_slices) ; 
    init_matrix_zero(sliceA, no_of_slices) ; 
    print_matrix(sliceA, no_of_slices, fp_matrix) ; 
    fprintf(fp_matrix, "# Initialize : matrix [%d] x [%d] |sliceB|  with Zero ... \n", 1, no_of_slices) ; 
    init_matrix_zero(sliceC, no_of_slices) ; 
    print_matrix(sliceC, no_of_slices, fp_matrix) ; 
    fprintf(fp_matrix, "# Initialize : matrix [%d] x [%d] |allC|  with Zero ... \n", nx, ny) ; 
    init_matrix_zero(allC, nx*ny) ; 
    print_matrix(allC, nx*ny, fp_matrix) ; 
     
//  MANUAL execution : calculate |allC| and infinity norm for resulting |C|
    fprintf(stdout, "# RESULTS : manual Straight-forward IJK calculation ... \n") ; 
    fprintf(fp_matrix, "\n# RESULTS : manual Straight-forward IJK calculation ... \n") ; 
    fprintf(fp_matrix, "# Initialize matrix [%d] x [%d] |C| for Straight-forward IJK BLAS/ATLAS computation ... \n", nx, ny) ; 
    init_matrix_zero(allC, nx*ny) ; 
    print_matrix(allC, nx*ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double manual_inf_norm = multipy_abc_manual(allA, allB, allC, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double manual_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6 ; 
    fprintf(fp_matrix, "# |C| : [%d] x [%d] matrix computed values : MANUAL complex ... \n", nx, ny) ; 
    print_matrix(allC, nx*ny, fp_matrix) ; 
    fprintf(fp_matrix, "# |C| : matrix calculated in %f seconds and has infinity norm of %g  ... \n", manual_elapsed, manual_inf_norm) ; 

//  CALCULATION :: |C| using DGEMM
    fprintf(stdout, "# RESULTS : Straight-forward IJK BLAS/ATLAS computation ...\n") ; 
    fprintf(fp_matrix, "\n# RESULTS : Straight-forward IJK BLAS/ATLAS computation ...\n") ; 
    fprintf(fp_matrix, "# Initialize matrix [%d] x [%d] |C| for Straight-forward IJK BLAS/ATLAS computation ... \n", nx, ny) ; 
    init_matrix_zero(allC, nx*ny) ; 
    print_matrix(allC, nx*ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    double dgemm_inf_norm = multipy_abc_cblas(allA, allB, allC, nx, ny) ; 
    gettimeofday(&tv2, &tz) ; 
    double dgemm_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6 ; 
    fprintf(fp_matrix, "\n# |C| : [%d] x [%d] matrix computed values using Straight-forward IJK BLAS/ATLAS computation ... \n", nx, ny) ; 
    print_matrix(allC, nx*ny, fp_matrix) ; 
    fprintf(fp_matrix, "# |C| : calculated in %f seconds and has infinity norm of %g  ... \n", dgemm_elapsed, dgemm_inf_norm) ; 
    

//  PTHREAD execution : Commence pthread calculation
    fprintf(stdout, "RESULTS : pthread computation ... \n") ; 
    fprintf(fp_matrix, "# RESULTS : pthread computation ... \n") ; 
    fprintf(fp_matrix, "# Initialize matrix [%d] x [%d] |C| for Straight-forward IJK BLAS/ATLAS computation ... \n", nx, ny) ; 
    init_matrix_zero(allC, nx*ny) ; 
    print_matrix(allC, nx*ny, fp_matrix) ; 
    gettimeofday(&tv1, &tz) ; 
    pthread_t *working_thread = malloc(np * (sizeof(pthread_t) )) ; 
    pthread_mutex_t *mutex_dot_product = malloc((sizeof(pthread_mutex_t))) ; 
    norm_sum_t *thread_norm_sum_data = malloc(np * (sizeof(norm_sum_t ))) ; 
    double *norm_sum_matrix = malloc( no_of_slices * (sizeof(double))) ; 
    void *status ; 
    double pthread_norm =0.0 ; 
    
    pthread_mutex_init(mutex_dot_product, NULL) ; 
    fprintf(stdout, "# DEBUG : \tCalculate infinity norm of |allB| and |allA| using |sliceA| writing to |sliceC|...\n") ; 
    gettimeofday(&tv1, &tz) ; 
    for(i=0 ; i<np ; i++) 
    { 
        fprintf(stdout, "# DEBUG :\tCreate working thread [%d] and send to norm_row_calculation function ...\n", i) ; 
        thread_norm_sum_data[i].my_sliceA = allA + (i * no_of_slices) ; 
        thread_norm_sum_data[i].my_sliceC =  allC + (i * no_of_slices) ; 
        thread_norm_sum_data[i].my_allB = allB ; 
        thread_norm_sum_data[i].my_norm_summation = norm_sum_matrix ; 
        thread_norm_sum_data[i].my_norm_sumC = -1.0 ; 
        thread_norm_sum_data[i].mutex = mutex_dot_product ; 
        thread_norm_sum_data[i].my_no_of_slices = no_of_slices ; 
        thread_norm_sum_data[i].my_nx = nx ; 
        thread_norm_sum_data[i].my_np = np ;       
        // double *my_sliceA ; 
        // double *my_sliceC ; 
        // double *my_allB ; 
        // double *my_norm_summation ; 
        // double my_norm_sumC ; 
        // pthread_mutex_t *mutex ; 
        // int my_no_of_slices ; 
        // int my_nx ; 
        // int my_np ; 
        pthread_create(&working_thread[i], NULL, norm_row_calculation, (void*)&thread_norm_sum_data[i]) ; 
    
        double pthread_row_norm = thread_norm_sum_data[i].my_norm_sumC ; 
        fprintf(stdout, "# DEBUG :\tbefore if - pthread_row_norm = [%g] \n", pthread_row_norm ) ; 
        fprintf(stdout, "# DEBUG :\tbefore if - pthread_row_norm = [%g]\n", pthread_row_norm) ;    
        pthread_norm = ( pthread_norm < pthread_row_norm ) ? pthread_row_norm : pthread_norm ; 
        fprintf(stdout, "# DEBUG :\tafter if - pthread_row_norm = [%g]\n", pthread_row_norm ) ;    
        fprintf(stdout, "# DEBUG :\tafter if - pthread_row_norm = [%g]\n", pthread_row_norm) ;    
    } 
     
    for(i=0 ; i<np ; i++) 
    {
        fprintf(stdout, "# DEBUG :\tImplement pthread_join for working thread[%d] ... \n", i) ; 
        pthread_join(working_thread[i], &status) ; 
    }
    printf("\tCalculate infintity norm of matrix C containing dot product of A and B using <mutex>...\n") ; 
      
//  program execution time taken : end
     gettimeofday(&tv2, &tz) ; 
     double pthread_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6 ; 

//  OUTPUT :: results to stdout & .dat file : matrix size || no of processors ||infinity norm || Time/manual || infinity norm || Time / dgemm
    fprintf(stdout, "# \t\t|Matrix|  |Threads|  Time/manual  Inf Norm/manual  Time/dgemm   Inf Norm/dgemm   Time/pthreads  Inf Norm/pthreads\n") ; 
    fprintf(stdout, "# Results: \t%d \t%d \t%f \t%g \t%f \t%g \t%f \t%g \n", nx, np, manual_elapsed, manual_inf_norm, dgemm_elapsed, dgemm_inf_norm, pthread_elapsed, pthread_norm) ; 
    fprintf(fp_timing, "%d \t%d \t%f \t%g \t%f \t%g \t%f \t%g \n", nx, np, manual_elapsed, manual_inf_norm, dgemm_elapsed, dgemm_inf_norm, pthread_elapsed,  pthread_norm) ; 

//  CLEANUP & close files
    fprintf(stdout, "# CLEAN-UP ... \n") ; 
    deallocate_matrix_memory(allA) ; 
    deallocate_matrix_memory(allB) ; 
    deallocate_matrix_memory(allC) ;    
    deallocate_matrix_memory(sliceA) ; 
    deallocate_matrix_memory(sliceC) ; 
     free(working_thread) ; 
     free(thread_norm_sum_data) ; 
     pthread_mutex_destroy(mutex_dot_product) ; 
     free(mutex_dot_product) ; 
    fclose(fp_matrix) ; 
    fclose(fp_timing) ; 

  return 0 ; 
}
