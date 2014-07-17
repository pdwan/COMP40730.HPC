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
#include <math.h>
#include <string.h>
#include <mpi.h>

int nx, ny, np, rank;

void usage () 
{
    fprintf(stdout,"\nUSAGE :\t<program name> [<-r>|<-i>] [N] <matrix contents file>.txt <timing file>.dat \n");
    fprintf(stdout,"\nTO :\tCalculate |allC| = |allA| x |allB| using MPI and also calculate infinity norm of |allC|. \n");    
    fprintf(stdout,"\nWHERE :\t1. \t<-r> \tinitialize |allA| & |allB| with _random_ numbers and |allC| with '0' \n");
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
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program : \tA4-mpi-1D \n# where : \t.dat contains timing data & .txt contains matrix values \n# \n");
    fprintf(fp, "# Summary of values added to each matrix - retained for later reference and validation \n# \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------  \n ");
}

void initialize_data_file_contents (FILE *fp) 
{
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# Program : \tA4-mpi-1D \n# where : \t.dat contains timing data & .txt contains matrix values \n");
    fprintf(fp, "#  --------------------------------------------------------------------------------------------------\n# \n");
    fprintf(fp, "# |Matrix| \t|Processors| \tTime/mpi \tInf Norm/mpi\n# \n");
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of main

int main (int argc, char *argv[])
{

//  define variables
    int increment_or_random = 0;    // increment enabled by default
    int max_num_args = 5;
    char filename_matrix[50];
    char filename_timing[50];
    int MAXN = 100;
    double sum, start;
    double *segmentA, *segmentB, *segmentC, *allB, *allC ;  
    double mpi_norm, mpi_row_norm,  mpi_elapsed; 
    int ni, nj, nk, my_nx;

//  MPI Initialization
    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          
//  CLI PARAMETERS :: validate and initialize
    if (argc != max_num_args) 
    {
        fprintf(stderr, "\nERROR: \t<number of arguments> [%d] : is invalid, less than <default> [%d].\n", argc, max_num_args);      
        MPI_Finalize();  
        usage();
    }  
//  random or increment initialization of matrices 
    char init_type[3];
    strncpy(init_type, argv[1], 2);
    init_type[3] = '\0';
    if ( strncmp( init_type , "-i" , 2 ) == 0 ) 
    {   
        increment_or_random = 0; 
    } else if ( strncmp( init_type , "-r" ,  2 ) == 0 ) 
    {   
        increment_or_random = 1; 
    } else
    { 
        fprintf(stderr, "\nERROR :\tInvalid entry : [%.2s] for random : '-r' or incremental : '-i'.\n", init_type);        
        MPI_Finalize();
        usage();
    }
//  matrix size
    nx = atoi(argv[2]);                                     
    if (nx > MAXN)
    {
        fprintf(stderr, "\nWARNING: \tMatrix size entered <nx> [%d]  too large, now set to [%d]. \n", nx, MAXN); 
        nx = MAXN;
    }    
    ny = nx;        
//  validate number of processors 
    if ((nx % np) != 0)
    {
        fprintf(stderr, "\nWARNING : \t<np> [%d] : number of processors must divide evenly into <nx> [%d] : matrix size. Existing. \n", np, nx) ; 
        MPI_Finalize();  
        usage(); 
    }
    if (np > nx)
    {
        fprintf(stderr, "\nWARNING: \t<np> [%d] : number of processors must be less than <nx> [%d] : matrix size. Exiting. \n", np, nx) ; 
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
    fprintf(fp_matrix, "# \n# RUNNING : \t%s %.2s %d %d \n", argv[0], init_type, nx, np);
    fprintf(stdout, "\n# RUNNING : \t%s %.2s %d %d \n", argv[0], init_type, nx, np);

//  ALLOCATE & INITIALIZE :: matrices and output results to matrix file for reference
    my_nx = nx / np;

    fprintf(stdout,"# ALLOCATE :\t|segmentA|, |segmentB|, |segmentC| and |allB| ... \n");
    fprintf(fp_matrix, "# ALLOCATE :\t|segmentA|, |segmentB|, |segmentC| and |allB| ... \n");
    segmentA = malloc(my_nx*ny*sizeof(double));
    segmentB = malloc(my_nx*ny*sizeof(double));
    segmentC = malloc(my_nx*ny*sizeof(double));
    allB = malloc(nx*ny*sizeof(double));
    fprintf(stdout,"# INITIALIZE :\tmatrices ... \n");
    fprintf(fp_matrix,"# INITIALIZE :\tmatrices ... \n");
    if (increment_or_random == 1) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |segmentA| using random mod 10 ... \n", my_nx, ny) ; 
        initialize_matrix_random(segmentA, my_nx, ny);
        print_matrix(segmentA, my_nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t<%d> x <%d> matrix |segmentB| using random mod 10 ... \n", my_nx, ny) ; 
        initialize_matrix_random(segmentB, my_nx, ny);
        print_matrix(segmentB, my_nx, ny, fp_matrix);
    } else if (increment_or_random == 0) 
    {
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |segmentA| using incremental <column> value + 1... \n", my_nx, ny) ; 
        initialize_matrix_increment(segmentA, my_nx, ny);
        print_matrix(segmentA, my_nx, ny, fp_matrix);
        fprintf(fp_matrix, "# \t\t<%d> x <%d> matrix |segmentB| using incremental <column> value + 1 ... \n", my_nx, ny) ; 
        initialize_matrix_increment(segmentB, my_nx, ny);
        print_matrix(segmentB, my_nx, ny, fp_matrix);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) 
    {
        start = MPI_Wtime(); 
    }
    for (ni=0; ni<np; ni++) 
    {
        MPI_Gather(segmentB, my_nx*ny, MPI_DOUBLE, allB, my_nx*ny, MPI_DOUBLE, ni, MPI_COMM_WORLD);
    }
    for (ni=0; ni<my_nx; ni++) 
    {
        for (nj=0; nj<ny; nj++) 
        {
            double sum = 0.0;
            for (nk=0; nk<nx; nk++)
            {
                sum+= (segmentA[(ni*nx)+nk]) * (allB[(nk*nx)+nj]);
//                fprintf(stdout,"DEBUG : \t = segmentA[(ni*nx)+nk] * allB[(nk*nx)+nj] : sum = segmentA[%d] * allB[%d] : %g += %g * %g \n",  (ni*nx)+nk, (nk*nx)+nj,sum , segmentA[(ni*nx)+nk], allB[(nk*nx)+nj] );
            }
            segmentC[(ni*nx)+nj] = sum;
//            fprintf(stdout,"DEBUG : \tsegmentC[(ni*nx) + nj] = sum : segmentC[%d] = sum : %g = %g \n", (ni*nx) + nj, segmentC[(ni*nx) + nj], sum);
         }
    }   
    fprintf(fp_matrix, "# VALIDATE :\t<%d> x <%d> matrix |allB| contents ... \n", nx, ny) ; 
    print_matrix(allB, nx, ny, fp_matrix);
    free(allB);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)  
    {
        fprintf(fp_matrix, "# ALLOCATE :\t<%d> x <%d> matrix |allC|  ... \n", nx, ny) ; 
        allC = malloc(nx*ny*sizeof(double));
    }  
    MPI_Gather (segmentC, my_nx*ny, MPI_DOUBLE, allC, my_nx*ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
        mpi_elapsed = MPI_Wtime() - start; 
        fprintf(stdout, "# RESULTS : \tMPI computation ...\n") ; 
        fprintf(stdout, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n", mpi_elapsed, mpi_norm) ;     
        fprintf(fp_matrix, "# RESULTS : \tMPI computation ...\n") ; 
        fprintf(fp_matrix, "# \t\tComputed Matrix [%d] x [%d] |allC| ... \n", nx, ny) ; 
        print_matrix(allC, nx, ny, fp_matrix) ; 
        fprintf(fp_matrix, "# \t\tMatrix |allC| calculated in [%f] seconds and has infinity norm of [%g] ... \n\n", mpi_elapsed, mpi_norm) ; 

//  OUTPUT :: results to stdout & .dat file : |Matrix| || |Processors| || Time / mpi || Inf norm/mpi
        fprintf(stdout,"# SUMMARY:\t|Matrix|  |Processors| Time/mpi Inf Norm/mpi\n");
        fprintf(stdout,"\t\t%d \t%d \t%f \t%g \n", nx, np, mpi_elapsed, mpi_norm);
        fprintf(fp_timing, "%d \t%d \t%f \t%g \n", nx, np, mpi_elapsed, mpi_norm);
    }

//  CLEANUP & close files
    fprintf(stdout,"# CLEAN-UP ... \n");
    if (rank == 0) 
    {
        free(allC);    
    }
    free(segmentA);
    free(segmentB);
    free(segmentC);    
    MPI_Finalize();  
    fclose(fp_matrix);
    fclose(fp_timing);
    return 0;
}
