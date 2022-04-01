#ifndef _ALLVEC_
#define _ALLVEC_

double **dmatrix(int nr, int nc);
int **imatrix(int nr, int nc);
void free_dmatrix(double **dmtx, int M);
void free_imatrix(int **imtx, int M);

#endif


 
