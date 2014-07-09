/***********************************************************************************************
  			ProMCore 2008 
		    February 05 - 09, 2008


 Example 2.3	: pthread-infinitynorm-rowwise.c
	

 Objective  	: Inifinity Norm of Matrix using Row Wise splitting.

 Input      	: Class Size,
         
                  Number of threads.
            
 Output     	: Infinity Norm of Matrix (Using Row-Wise splitting).

 Created	: January 2008                                              
               	  National PARAM Supercomputing Facility,C-DAC,Pune.        

 E-mail      	: betatest@cdac.in  

***********************************************************************************************/

#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include<stdlib.h>
#include<sys/time.h>

#define MAXTHREADS 8

/* pthread mutex and condition variable declaration */
pthread_mutex_t  mutex_norm = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  CurRow_norm = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_Row = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t count_threshold_cv = PTHREAD_COND_INITIALIZER;

/*global variable declaration */
pthread_mutex_t  *mutex_Res;

double   NormRow = 0 ;
double   *Res;
int      dist_row,dist_col;
int      iteration;

int      rmajor,cmajor;
int      perfect_square;
int      row, col, currentRow = 0;
double   **InMat ;
int      numberOfThreads;
double   infinitynorm_row;

/* Thread callback function */
void  * doRowWise(int myId)
{
       int   CurRow = 0;      
       int   iCol,myRowSum;
       int   mynorm;

       for (CurRow = ((myId - 1) * dist_row); CurRow <= ((myId * dist_row) - 1); CurRow++)
       {
         
            myRowSum = 0;
            for(iCol = 0 ;iCol < col; iCol++)
            myRowSum += InMat[CurRow][iCol];

            if(mynorm < myRowSum )
            mynorm = myRowSum;
       }

       pthread_mutex_lock(&mutex_norm);
       {
             if (NormRow < mynorm)
              NormRow = mynorm;
       }
       pthread_mutex_unlock(&mutex_norm);

       pthread_exit(NULL);

}

/* Main function starts*/
main(int argc, char *argv[])
{
/* variable declaration */
  int            i, j,p,q,CLASS_SIZE,THREADS;
  int            first_value_row,diff,matrix_size,ret_count;
  double         time_start, time_end,memoryused=0.0;
  struct         timeval tv;
  struct         timezone tz;
  int            counter, irow, icol, numberOfThreads,ithread;  
  FILE           *fp;
  char            * CLASS;

  pthread_t      *threads_row;
  pthread_attr_t ptr;

  
  printf("\n\t\t---------------------------------------------------------------------------");
  printf("\n\t\t Centre for Development of Advanced Computing (C-DAC):  Febraury-2008");
  printf("\n\t\t Email : betatest@cdac.in");
  printf("\n\t\t---------------------------------------------------------------------------");
  printf("\n\t\t Objective : Dense Matrix Computations (Floating-Point Operations)\n ");
  printf("\n\t\t Computation of Infinity Norm of a Square Matrix - Rowwise partitioning;");
  printf("\n\t\t..........................................................................\n");
   /*Argument validation check*/
  if( argc != 3 ){

     printf("\t\t Very Few Arguments\n ");
     printf("\t\t Syntax : exec <Class-Size> <Threads>\n");
           exit(-1);
   }
   else {
      CLASS = argv[1];
      THREADS = atoi(argv[2]);
   }
   if( strcmp(CLASS, "A" )==0){
       CLASS_SIZE = 1024;
   }
   if( strcmp(CLASS, "B" )==0){
       CLASS_SIZE = 2048;
   }
   if( strcmp(CLASS, "C" )==0){
         CLASS_SIZE = 4096;
   }

  numberOfThreads = THREADS;
   matrix_size = CLASS_SIZE;
   row = matrix_size;
   col = matrix_size;

   printf("\n\t\t Input Parameters :");
   printf("\n\t\t CLASS           :  %c",*CLASS);
   printf("\n\t\t Matrix Size     :  %d",matrix_size);
   printf("\n\t\t Threads         :  %d",numberOfThreads);
   printf("\n");

  if ((numberOfThreads != 1) && (numberOfThreads != 2) && (numberOfThreads != 4) && (numberOfThreads != 8))
   {
       printf("\n Number of Threads must be 1 or 2 or 4 or 8. Aborting ...\n");
       exit(0);
   }

  if(numberOfThreads > row)
   {
       printf("\nNumber of threads should be <= %d",row);
       exit(0);
   }

  /*....Memory Allocation....*/
  
  InMat = (double **) malloc(sizeof(double) * row);
  for (i = 0; i < row; i++)
  InMat[i] = (double *) malloc(sizeof(double) * col);
  if(InMat == NULL)
  {
   	printf("\n Not sufficient memory to accomadate InMat1");
   	exit(0);
  } 
  memoryused = (row * col * sizeof(double));


  /* Matrix  and Vector  Initialization */
   for (i=0; i<row;i++)
   {
   	for (j = 0; j<col;j++)
   	{
   		InMat[i][j] = (double)(i + j);
   	}
   }

  
 /* InfinityNorm Of InMat1 Matrix */

   /* Row Wise Partitioning */

  dist_row = row/numberOfThreads;

  threads_row = (pthread_t *) malloc(sizeof(pthread_t) * numberOfThreads);

   /*Taking start time*/
  gettimeofday(&tv, &tz);
  time_start = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;

  /*Thread creation */
  pthread_attr_init(&ptr);
  for (ithread = 0; ithread < numberOfThreads; ithread++)
  {
  	ret_count=pthread_create(&threads_row[ithread], &ptr, (void *(*) (void *)) doRowWise, (void *) (ithread+1));
	if(ret_count)
    	{
		printf("\n ERROR : Return code from pthread_create() is %d ",ret_count);
		exit(-1);
    	}
  }
  /*Joining threads */
  for (ithread = 0; ithread < numberOfThreads; ithread++)
  {
  	ret_count=pthread_join(threads_row[ithread], NULL);
 	if(ret_count)
    	{
		printf("\n ERROR : Return code from pthread_join() is %d ",ret_count);
		exit(-1);
    	}
  }

  printf("\n\t\t Row Wise partitioning - Infinity Norm of a Square Matrix.....Done ");

  pthread_attr_destroy(&ptr);


   free(InMat);
   free(threads_row);
   gettimeofday(&tv, &tz);
   time_end = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
   printf("\n");
   printf("\n\t\t Memory Utilized       :  %lf MB",(memoryused/(1024*1024)));
   printf("\n\t\t Time in  Seconds (T)  :  %lf",(time_end - time_start));
   printf("\n\t\t..........................................................................\n");


}

