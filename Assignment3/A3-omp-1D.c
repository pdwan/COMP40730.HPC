/* 
* Paula Dwan
* 13208660 : COMP40730 : Assignment 3
*
* Open MP program computing norm of product of two nxn matrices
* 1/     Granularity of the product : 
*           (ii)      One step algorithim. No intermediate resulting matrix
* 2/     Partitioning Scheme :
*           (a)      Left matrix is horizonitally partitioned
* 3/     Matrix norm to be computed :
*           (b)      Maximum absolute row sum norm (aka infinity-norm)
*                     calculate max value of sum of each row (or product of same) 
*                     ||A|| infintity = max SUM | Aij |
*                     where   max     :==> 0 <= i <=n
*                                  SUM   :==>  j =0 to n-1      
*
*/          

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <cblas.h>
#include <math.h>
#include <omp.h>

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
// function to obtain value from user for matrix size, number of processes.

  return 0;
}

