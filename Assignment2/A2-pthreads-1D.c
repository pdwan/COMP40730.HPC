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

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <cblas.h>
#include <math.h>
#include <string.h>

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
typedef struct 
{ 
     double *my_sliceA; 
     double *my_sliceC; 
     double *my_allB; 
     double *my_norm_summation; 
     double my_norm_sumC; 
     pthread_mutex_t *mutex; 
     int my_ns_size;
     int my_nx; 
     int my_np;
} norm_summation_t;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void usage () {
    fprintf(stdout,"\nUSAGE : \t<program name> [<-r>|<-i>] [N] [P] <matrix contents file>.txt <timing file>.dat \n");
    fprintf(stdout,"\nTO : \t\tCalculate |C| = |A| x |B| using pthreads and also calculate infinity norm of |C| \n");    
    fprintf(stdout,"\nWHERE :");
    fprintf(stdout,"\t1. \t<-r>\tinitialize |A| & |B| with _random_ numbers and |C| with '0' \n");
    fprintf(stdout,"\t   \t<-i>\tinitialize |A| & |B| _incrementally_ with <column> value and |C| with '0' \n");
    fprintf(stdout,"\t2. \t[N] \tmatrix size, if invalid defaults to [${maxMatrixSize}] \n");
    fprintf(stdout,"\t3. \t[P] \tnumber of threads, i.e.: processors when (i) [P] less than [N] and (ii) [N] mod [P] = 0 \n\t\tIf invaldi, set to [${max}]");
    fprintf(stdout,"\t4. \t<matrix contents file>.txt\n\t\tname of .txt file to store values of matrices |A| |B| & |C| \n");
    fprintf(stdout,"\t5. \t<timing .dat file> .dat \n\t\tname of timing data file to containing calculation time for each iteration \n\n");
    exit(0);
}

double* allocate_memory_matrix (int rows_cols)
{
    double* l_matrix;  
    l_matrix = (double *) malloc(rows_cols*sizeof (double));  
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

void print_matrix(double *matrix1d, int rows_cols, FILE *fp)
{
    int ni;
    int count =1;
    for (ni=0; ni<(rows_cols); ni++)
    {
        fprintf(fp, "%g\t", matrix1d[ni]); 
        if (count == sqrt(rows_cols) ) 
        {
            fprintf(fp, "\n");
            count =1;
        } else 
        {
            count ++;    
        }
    }
}

void init_matrix_random(double *matrix1d, int rows_cols)
{
    int ni, mod_rand=10;
    for (ni=0; ni<(rows_cols); ni++)
    {
        matrix1d[ni]= (double) ((rand()%mod_rand) + 1);
    }
}

void init_matrix_increment(double *matrix1d, int rows_cols)
{
    int ni;    
    int n_rows_cols= sqrt(rows_cols) ;
    for (ni=0; ni<(rows_cols); ni++)
    {
        matrix1d[ni]= ( ( ni % n_rows_cols ) + 1);
    }
}

void init_matrix_zero(double *matrix1d, int rows_cols)
{
    int ni;
    for (ni=0; ni<(rows_cols); ni++)
    {
        matrix1d[ni]=  0.0;
    }
}

int validate_if_file_exists(char * fn)
{
    FILE* fp = fopen(fn, "r");
    if (fp != NULL)
    {
        fclose(fp);
        return 1;   // file exists
    }
    return 0;       // file does not exist
}

void init_matrix_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA2-pthreads-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n");
}

void init_timing_file_contents(FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program :\tA2-pthreads-1D \n# where :\t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n");
    fprintf(fp, "# |Matrix| \t|Threads| \tTime/manual \tInf Norm/manual \tTime/dgemm \tInf Norm/dgemm \n# \n");
}

//  function to calculate infinity norm and provide values for |C| : manual
double multipy_ABC_complex(double *matrix_a, double *matrix_b, double *matrix_c, int rows, int cols, FILE *fp)
{
    int ni, nj, nk;
    
    for (ni=0; ni<rows; ni++)
    {
        for (nj=0; nj<cols; nj++)
        {
            double sum = 0.0;
            for (nk=0; nk<rows; nk++)
            {
                sum+= (matrix_a[(ni*rows)+nk]) * (matrix_b[(nk*rows)+nj]);
            }
            matrix_c[(ni*rows)+nj] = sum;
        }
    }
double manual_row_norm =0.0, manual_inf_norm =0.0;
    for (ni=0; ni<rows; ni++)
    {
        manual_row_norm =0.0;
        for (nj=0; nj<rows; nj++)
        {
            manual_row_norm += matrix_c[(ni*rows) +nj];
        }
        manual_inf_norm = ( manual_inf_norm < manual_row_norm ) ? manual_row_norm : manual_inf_norm;
    }
    return manual_inf_norm;
}

//  function to calculate infinity norm and provide values for |C|
void *norm_summation_calculation(void *args)
{
     norm_summation_t *norm_summation_data; 
     int i;  
     int ALPHA=1.0;
     int BETA=0.0;
     
     fprintf(stdout, "# Calculate infinity norm of |allB| and |allA| using |sliceA| writing to |sliceC|...\n");
     norm_summation_data = args; 
     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, norm_summation_data->my_nx, norm_summation_data->my_nx, norm_summation_data->my_nx, ALPHA, norm_summation_data->my_sliceA, norm_summation_data->my_nx, norm_summation_data->my_allB, norm_summation_data->my_nx, BETA, norm_summation_data->my_sliceC, norm_summation_data->my_nx);
     fprintf(stdout, "# Evaluate sum of local |sliceC| using <mutex>...\n");
     pthread_mutex_lock(norm_summation_data->mutex);
     fprintf(stdout, "# ................... > norm_summation_data->my_ns_size = %d \n ", (int)norm_summation_data->my_ns_size);
     fprintf(stdout, "# value \t\t  max");
     int iteration = (int) (norm_summation_data->my_ns_size / norm_summation_data->my_nx);
     int counter = 0;
     for(i=0; i<iteration; i++) 
     {
          norm_summation_data->my_norm_summation[counter] += norm_summation_data->my_sliceC[i];
          if ((counter % (int) norm_summation_data->my_nx) > 0)
          { 
               counter ++;
          } else 
          {
               counter =0;
          }
          if (norm_summation_data->my_norm_sumC < norm_summation_data->my_norm_summation[i]) 
          {
               norm_summation_data->my_norm_sumC = norm_summation_data->my_norm_summation[i] ;
          }
     fprintf(stdout, "# %g \t\t %g \n",  (double)norm_summation_data->my_norm_summation[i], norm_summation_data->my_norm_sumC);
     }  
     fprintf(stdout, "# Max absolute value in local |sliceC| is : %g\n", norm_summation_data->my_norm_sumC);     
     pthread_mutex_unlock(norm_summation_data->mutex); 
     pthread_exit(NULL); 
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// start of main
int main ( int argc, char *argv[] )
{
//  define variables
    struct timeval tv1, tv2;
    struct timezone tz;
    const int MAXPTHREADS = 100;
    const int MAXN = 1000;
    int i = 0;
    int increment_or_random = 5;
    int max_num_args = 6;
    char filename_matrix[50];
    char filename_timing[50];
              
//  start of program

    if ( argc != max_num_args )
    {
        usage();
    } 
    
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
    int np = atoi(argv[3]);             
    if ( (nx % np) != 0)
    {
        fprintf(stderr, "\nWARNING : \t<np> %d : number of threads must divide evenly into <nx> %d : matrix size. Using detaults. \n", np, nx);
        nx=MAXN;
        np= MAXPTHREADS;        
    }
    if (np > nx)
    {
        fprintf(stderr, "\nWARNING: \t<np> %d : number of threads must be less <nx> %d : matrix size. Using detaults. \n", np, nx);
        nx=MAXN;
        np= MAXPTHREADS;
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
        fp_matrix = fopen(filename_matrix, "a" );
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
        fp_timing = fopen(filename_timing, "a" );
    } 
    fprintf(fp_matrix, "RUNNING :\t%s %s %d %d %s %s \n", argv[0],cell_value_type,nx,np,filename_matrix,filename_timing );
    fprintf(stdout, "RUNNING : \t%s %s %d %d %s %s \n", argv[0],cell_value_type,nx,np,filename_matrix,filename_timing );
 
//  Calculate slice information using number of processors <np> and matrix size <nx>
     int ns_size = (nx*np);
     int ns_no = (nx / np); 
     fprintf(fp_matrix, "\tNo of slices = <matrix rows> x <number of processes> = [%d] x [%d] = [%d] ...\n", nx, np, ns_size);
     fprintf(fp_matrix, "\tSlice size       = <matrix rows> ÷ <number of processes> = [%d] ÷ [%d] = [%d] ...\n", nx, np, ns_no);
     fprintf(stdout, "\tNo of slices = <matrix rows> x <number of processes> = [%d] x [%d] = [%d] ...\n", nx, np, ns_size);
     fprintf(stdout, "\tSlice size       = <matrix rows> ÷ <number of processes> = [%d] ÷ [%d] = [%d] ...\n", nx, np, ns_no);
 
// Allocate memory for matrices allA, allB allC, sliceA and sliceC
     fprintf(stdout, "CREATE MATRICES |allA|, |allB|, |allC|, |sliceA| & |sliceC| ... \n");
     double *allA = allocate_memory_matrix(nx*ny);
     double *allB = allocate_memory_matrix(nx*ny);
     double *allC = allocate_memory_matrix(nx*ny);
     double *sliceA = allocate_memory_matrix(ns_size);
     double *sliceC = allocate_memory_matrix(ns_size);
 
 //  Initialize matrices A & B & C
    fprintf(stdout, "INITIALIZE MATRICES |allA|, |allB|, |allC|, |sliceA|  & |sliceC| ... \n");
    if (increment_or_random == 0) {
        fprintf(fp_matrix, "Initialize : matrix [%d] x [%d] |allA| using <column> values ... \n", nx, ny);
        init_matrix_random(allA, nx*ny);
        fprintf(fp_matrix, "Initialize : matrix [%d] x [%d] |allB| using <column> values ... \n", nx, ny);
        init_matrix_random(allB, nx*ny);
    } else if (increment_or_random == 1) 
    {
        fprintf(fp_matrix, "Initialize : matrix [%d] x [%d] |allA| using random numbers 1 to 10 ... \n", nx, ny);
        init_matrix_increment(allA, nx*ny);
        fprintf(fp_matrix, "Initialize : matrix [%d] x [%d] |allB|  using random numbers 1 to 10 ... \n", nx, ny);
        init_matrix_increment(allB, nx*ny);
    }
    print_matrix(allA, nx*ny,fp_matrix);
    print_matrix(allB, nx*ny,fp_matrix);
    fprintf(fp_matrix, "Initialize : matrix [%d] x [%d] |sliceA| with Zero ... \n", 1, ns_size);
    init_matrix_zero(sliceA, ns_size);
    print_matrix(sliceA, ns_size, fp_matrix);
    fprintf(fp_matrix, "Initialize : matrix [%d] x [%d] <sliceB>  with Zero ... \n",1, ns_size);
    init_matrix_zero(sliceC, ns_size);
    print_matrix(sliceC, ns_size, fp_matrix);
    fprintf(fp_matrix, "Initialize : matrix [%d] x [%d] |allC|  with Zero ... \n", nx, ny);
    init_matrix_zero(allC, nx*ny);
    print_matrix(allC, nx*ny, fp_matrix);
     
//  MANUAL execution : calculate |allC| and infinity norm for resulting |C|
    fprintf(stdout,"# RESULTS : complex manual calculation ... \n");
    gettimeofday(&tv1, &tz);
    double complex_inf_norm = multipy_ABC_complex(allA, allB, allC, nx, ny, fp_matrix);
    gettimeofday(&tv2, &tz);
    double complex_elapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
    fprintf(fp_matrix, "# |C| : <%d> x <%d> matrix computed values : MANUAL complex ... \n", nx, ny);
    print_matrix(allC, nx*ny, fp_matrix);
    fprintf(fp_matrix, "# |C| : matrix calculated in %f seconds and has infinity norm of %g  ... \n", complex_elapsed, complex_inf_norm);

// PTHREAD execution : Commence pthread calculation
    fprintf(stdout, "RESULTS : pthread computation ... \n");
    fprintf(fp_matrix, "# Initialize matrix <%d> x <%d> |C|, redone for pthread computation .. \n", nx, ny);
    init_matrix_zero(allC, nx*ny);
    print_matrix(allC, nx*ny, fp_matrix);       
    pthread_t *working_thread = malloc(np * (sizeof(pthread_t) ));
    pthread_mutex_t *mutex_dot_product = malloc(nx * (sizeof(pthread_mutex_t))); 
    norm_summation_t *thread_norm_summation_data = malloc(np * (sizeof(norm_summation_t))); 
    double norm_summation_matrix = malloc( nx * (sizeof(double)));
    void *status; 
    pthread_mutex_init(mutex_dot_product, NULL); 
    for(i=0; i<np; i++) 
    { 
        printf("\tCreate working thread[%d] and send to norm_summation_calculation function ...\n", i);
        thread_norm_summation_data[i].my_sliceA = allA + (i * ns_size);
        thread_norm_summation_data[i].my_sliceC =  allC + (i * ns_size);
        thread_norm_summation_data[i].my_allB = allB;
        thread_norm_summation_data[i].my_norm_summation = norm_summation_matrix;
        thread_norm_summation_data[i].my_norm_sumC = -1.0;
        thread_norm_summation_data[i].mutex = mutex_dot_product;
        thread_norm_summation_data[i].my_ns_size = ns_size;
        thread_norm_summation_data[i].my_nx = nx;
        thread_norm_summation_data[i].my_np = np;       
        double *my_sliceA; 
        double *my_sliceC; 
        double *my_allB; 
        double *my_norm_summation; 
        double my_norm_sumC; 
        pthread_mutex_t *mutex; 
        int my_ns_size;
        int my_nx; 
        int my_np;
        pthread_create(&working_thread[i], NULL, norm_summation_calculation,(void*)&thread_norm_summation_data[i]);
    } 
     
    for(i=0; i<np; i++) 
    {
        printf("\tImplement pthread_join for working thread[%d] ... \n", i);
        pthread_join(working_thread[i], &status); 
    }
    printf("\tCalculate infintity norm of matrix C containing dot product of A and B using <mutex>...\n");
      
// program execution time taken : end
     gettimeofday(&tv2, &tz);
     double pthread_elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;

//  OUTPUT :: results to stdout & .dat file : matrix size || no of processors ||infinity norm || Time/manual || infinity norm || Time / dgemm
    fprintf(stdout,"# \t\t|Matrix|  |Threads|  Time/manual  Inf Norm/manual  Time/pthreads  Inf Norm/pthreads\n");
    fprintf(stdout,"# Results: \t%d \t%d \t%fs \t%g \t%fs \t%g \n", nx, np, complex_elapsed, complex_inf_norm, pthread_elapsed,  norm_summation_data->my_norm_sumC);
    fprintf(fp_timing, "%d \t%d \t%fs \t%g \t%fs \t%g \n", nx, np, complex_elapsed, complex_inf_norm, pthread_elapsed,  norm_summation_data->my_norm_sumC);

//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    deallocate_matrix_memory(allA);
    deallocate_matrix_memory(allB);
    deallocate_matrix_memory(allC);    
    deallocate_matrix_memory(sliceA);
    deallocate_matrix_memory(sliceC);
     free(working_thread);
     free(thread_norm_summation_data);
     pthread_mutex_destroy(mutex_dot_product);
     free(mutex_dot_product);
    fclose(fp_matrix);
    fclose(fp_timing);

  return 0;
}

