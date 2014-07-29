#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>
#include <string.h>
#include <pthread.h>

double *A;
double *B;
double *C;
int n;
double matrix_norm;

typedef struct {
   double *b;
   double *c;
   int num_of_columns;
   pthread_mutex_t *mutex;
} matrix_slice;

void *matrix_slice_multiply(void *arg){
   matrix_slice *slice = arg;
   int i, j;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, slice->num_of_columns, n, 1.0, A, n, slice->b, n, 0.0, slice->c, n);

   // compute column norm of each slice
   double slice_norm = 0.0;
   for(j = 0; j < slice->num_of_columns; j++) {
      double column_sum=0.;
      for(i = 0; i < n; i++)
         column_sum += *(slice->c + i * n + j);

      if(column_sum>slice_norm)
         slice_norm=column_sum;
   }
   pthread_mutex_lock(slice->mutex);
   if (slice_norm>matrix_norm)
      matrix_norm=slice_norm;
   pthread_mutex_unlock(slice->mutex);

   pthread_exit(NULL);
}

int main(void) {
   int num_of_thrds, num_of_columns_per_slice;
   pthread_t *working_thread;
   matrix_slice *slice;
   pthread_mutex_t *mutex;
   int i = 0;

   printf ("Please enter matrix dimension n : ");
   scanf("%d", &n);

   printf ("Please enter number of threads : ");
   scanf("%d", &num_of_thrds);

   while (num_of_thrds > n) {
      printf("number of threads must not be greater than matrix dimension\n");
      printf ("Please enter number of threads : ");
      scanf("%d", &num_of_thrds);
   }
   // allocate memory for the matrices
   ///////////////////// Matrix A //////////////////////////
   A = (double *)malloc(n * n * sizeof(double));

   if (!A) {
      printf("memory failed \n");
      exit(1);
   }

   ///////////////////// Matrix B //////////////////////////
   B = (double *)malloc(n * n * sizeof(double));
   if (!B) {
      printf("memory failed \n");
      exit(1);
   }

   ///////////////////// Matrix C //////////////////////////
   C = (double *)malloc(n * n * sizeof(double));
   if (!C) {
      printf("memory failed \n");
      exit(1);
   }

   // initialize the matrices
   for (i = 0; i < n * n; i++) {
      A[i] = rand() % 15;
      B[i] = rand() % 10;
      C[i] = 0.;
   }

   clock_t t1 = clock();
   working_thread = malloc(num_of_thrds * sizeof(pthread_t));
   slice = malloc(num_of_thrds * sizeof(matrix_slice));
   mutex = malloc(sizeof(pthread_mutex_t));
   num_of_columns_per_slice = n / num_of_thrds;

   for(i = 0; i < num_of_thrds; i++){
      slice[i].b = B + i * num_of_columns_per_slice;
      slice[i].c = C + i * num_of_columns_per_slice;
      slice[i].mutex = mutex;
      slice[i].num_of_columns = (i == num_of_thrds - 1) ? n-i * num_of_columns_per_slice : num_of_columns_per_slice;
      pthread_create(&working_thread[i], NULL, matrix_slice_multiply, (void *)&slice[i]);
   }
   for(i = 0; i < num_of_thrds; i++)
      pthread_join(working_thread[i], NULL);

   clock_t t2=clock();
   printf("elapsed time: %f\n", (double)(t2 - t1)/CLOCKS_PER_SEC);

   printf("column sum norm is %f\n", matrix_norm);

   //deallocate memory
   free(A);
   free(B);
   free(C);
   free(working_thread);
   free(slice);

   return 0;
}
