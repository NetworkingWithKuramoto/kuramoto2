#include "type_handlers.h"


void nrerror(error_text)
char error_text[];
{
	/**
	 * @brief Prints an error message and exits the program.
	 *
	 * If a runtime error occurs (typically related to memory allocation), 
	 * this function provides a message detailing the error and terminates the program.
	 *
	 * @param error_text The error message to display.
	 */

	void exit();  

	printf("Numerical Recipes run-time error...\n");
	printf("%s\n",error_text);
	printf("...now exiting to system...\n");
	exit(1);
}

double **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	/**
	 * @brief Allocates a 2D matrix of doubles.
	 *
	 * Creates a matrix with rows ranging from `nrl` to `nrh` and columns ranging from `ncl` to `nch`.
	 *
	 * @param nrl The lower index of rows.
	 * @param nrh The upper index of rows.
	 * @param ncl The lower index of columns.
	 * @param nch The upper index of columns.
	 * @return A pointer to the allocated 2D matrix of type double.
	 */
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	/**
	 * @brief Allocates a 2D matrix of integers.
	 *
	 * Creates a matrix with rows ranging from `nrl` to `nrh` and columns ranging from `ncl` to `nch`.
	 *
	 * @param nrl The lower index of rows.
	 * @param nrh The upper index of rows.
	 * @param ncl The lower index of columns.
	 * @param nch The upper index of columns.
	 * @return A pointer to the allocated 2D matrix of type int.
	 */
	int i;
	int **m;

	m=(int  **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int  *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

double *vector(nl,nh)
int nl,nh;
{
	/**
	 * @brief Allocates a 1D vector of doubles.
	 *
	 * Creates a vector with indices ranging from `nl` to `nh`.
	 *
	 * @param nl The lower index of the vector.
	 * @param nh The upper index of the vector.
	 * @return A pointer to the allocated vector of type double.
	 */

	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
{
	/**
	 * @brief Allocates a 1D vector of integers.
	 *
	 * Creates a vector with indices ranging from `nl` to `nh`.
	 *
	 * @param nl The lower index of the vector.
	 * @param nh The upper index of the vector.
	 * @return A pointer to the allocated vector of type int.
	 */

	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

void free_matrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	/**
	 * @brief Frees the memory allocated for a 2D matrix of doubles.
	 *
	 * @param m Pointer to the matrix to be freed.
	 * @param nrl The lower index of rows.
	 * @param nrh The upper index of rows.
	 * @param ncl The lower index of columns.
	 * @param nch The upper index of columns.
	 */

	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
	/**
	 * @brief Frees the memory allocated for a 2D matrix of integers.
	 *
	 * @param m Pointer to the matrix to be freed.
	 * @param nrl The lower index of rows.
	 * @param nrh The upper index of rows.
	 * @param ncl The lower index of columns.
	 * @param nch The upper index of columns.
	 */

	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_vector(v,nl,nh)
double *v;
int nl,nh;
{
	/**
	 * @brief Frees the memory allocated for a 1D vector of doubles.
	 *
	 * @param v Pointer to the vector to be freed.
	 * @param nl The lower index of the vector.
	 * @param nh The upper index of the vector.
	 */

	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
	/**
	 * @brief Frees the memory allocated for a 1D vector of integers.
	 *
	 * @param v Pointer to the vector to be freed.
	 * @param nl The lower index of the vector.
	 * @param nh The upper index of the vector.
	 */

	free((char*) (v+nl));
}

