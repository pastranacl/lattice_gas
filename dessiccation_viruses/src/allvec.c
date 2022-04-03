/*****************************************************************************************/
/*                                      ALLVEC.C                                         */
/*                                                                                       */
/*  This library is analogous to nrutils from Numerical recipes in C. In summary, it     */
/*  contains the functions necessary for the allocation of 
/*****************************************************************************************/



#include <stdio.h>      
#include <stdlib.h>          


//**********************************************************************/
/*                               dmatrix                              */
/*  Allocate a double matrix with subscript range                     */                                                                 
/*  m[L rows][N cols]                                                 */
/*                                                                    */
/**********************************************************************/
double **dmatrix(int nr, int nc)
{
  double **m;

  // Allocate pointers to rows 
  m=(double **) malloc(nr*sizeof(double*));
  if (!m) 
  {
    fprintf(stderr, "Not enough memory for dmatrix, Err 1. \n");
    exit(1);
  }
  
  // Allocate rows and set pointers to them 
  for(int i=0; i<nr; i++)
        m[i] = (double *) malloc(nc*sizeof(double));

  // Return pointer to array of pointers to rows 
  return m;
}
 
 
/**********************************************************************/
/*                               imatrix                              */
/*  Allocate an int matrix with subscript range                       */                                                                 
/*  m[L rows][N cols]                                                 */
/*                                                                    */
/**********************************************************************/
int **imatrix(int nr, int nc)
{
  int **m;

  // Allocate pointers to rows 
  m=(int **) malloc(nr*sizeof(int*));
  if (!m) 
  {
    fprintf(stderr, "Not enough memory for imatrix, Err 1. \n");
    exit(1);
  }
  
  // Allocate rows and set pointers to them 
  for(int i=0; i<nr; i++)
        m[i] = (int *) malloc(nc*sizeof(int));

  // Return pointer to array of pointers to rows 
  return m;
}
 

 
/************************************************************************/
/*                             free_imatrix                             */
/* Free the memory for the matrix of M rows allocated with imatrix      */
/************************************************************************/
void free_imatrix(int **imtx, int M)
{
    for(int i=0; i<M;i++)
        free(imtx[i]);
    free(imtx); 
}


/************************************************************************/
/*                             free_dmatrix                             */
/* Free the memory for the matrix of M rows allocated with dmatrix      */
/************************************************************************/
void free_dmatrix(double **dmtx, int M)
{
    for(int i=0; i<M;i++)
        free(dmtx[i]);
    free(dmtx); 
}

