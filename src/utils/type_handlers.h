#ifndef TYPE_HANDLERS_H
#define TYPE_HANDLERS_H

#include <stdio.h>
#include <stdlib.h>

// Error handling function
void nrerror(char error_text[]);

// Matrix allocation functions
double **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);

// Vector allocation functions
double *vector(int nl, int nh);
int *ivector(int nl, int nh);

// Freeing memory functions
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_vector(double *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);

#endif // TYPE_HANDLERS_H
