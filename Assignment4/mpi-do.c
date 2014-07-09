/* MPI easy one */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

const int size = 1000;

float a [size][size];
float b [size][size];
float c [size][size];

void multiply (int istart, int iend) 
{
    int i=0; 
    int j=0; 
    int k=0;
    for (i =istart; i <=iend; ++i) 
    {
        for (j =0; j <size; ++j) 
        {
            for (k =0; k <size; ++k) 
            {
                c[i][j] += a[i][k] * b{k][j] ;
            }
        }
    }
}

int main ( int argc, char *argv[] )
{
    int rank, nproc;
    int istart, iend;
    int i=0; 
    int j=0; 

    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if ( rank == 0 ) 
    {
        for (i =istart; i <=iend; ++i) 
        {
            for (j =0; j <size; ++j) 
            {
                a[i][j] = (float) (i + j);
                b[i][j] = (float) (i - j);
                c[i][j] = 0.0f;
            }
        }
    }
    
    MPI_Bcast(a, size*size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, size*size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(c, size*size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    istart = (size / nproc) * rank;
    iend = (size / nproc) * (rank +1) - 1;
    
    multiply(istart, iend);
    
    MPI_Gather(c+ (size/nproc * rank), size/size/nproc, MPI_FLOAT,
        c+ (size/nproc * rank), size/size/nproc, MPI_FLOAT, 0, MPI_COMM_WORLD);
        
    if ( rank == 0) 
    {
        if ( (size % nproc) > 0 ) 
        {
        multiply( ((size/nproc)*nproc), (size-1) );
        }
    }        
    
    MPI_Finalize();
    exit 0;
}
