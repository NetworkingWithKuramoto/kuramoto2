#include "rk4.h"

void rk4(y,dydx,n,x,h,yout,derivs)
double y[],dydx[],x,h,yout[];
void (*derivs)();	/* ANSI: void (*derivs)(double,double *,double *); */
int n;
{
	/**
	* @brief Performs a single step of the fourth-order Runge-Kutta integration.
	* 
	* This function integrates a system of ordinary differential equations (ODEs) using 
	* the classic fourth-order Runge-Kutta (RK4) method over one time step of size `h`.
	* 
	* @param[in] y      An array of size `n+1` containing the initial values of the dependent variables.
	*                   The array index starts from 1, not 0.
	* @param[in] dydx   An array of size `n+1` containing the initial derivatives of `y` at `x`.
	*                   The array index starts from 1, not 0.
	* @param[in] n      The number of equations to be integrated.
	* @param[in] x      The independent variable (e.g., time) at the start of the step.
	* @param[in] h      The step size for integration.
	* @param[out] yout  An array of size `n+1` where the function stores the updated values 
	*                   of `y` after advancing by one step `h`.
	* @param[in] derivs A function pointer to the derivative function. 
	*                   The function `derivs` should have the following signature:
	*                   ```c
	*                   void derivs(double x, double y[], double dydx[]);
	*                   ```
	*                   It computes the derivatives of `y` at `x`.
	* 
	* @details
	* - The function uses temporary arrays `dym`, `dyt`, and `yt` to store intermediate steps 
	*   of the RK4 method.
	* - The RK4 method calculates a weighted average of four increments:
	*     - The slope at the beginning of the interval (`k1`).
	*     - The slope at the midpoint of the interval, using `k1` (`k2`).
	*     - The slope at the midpoint, but using `k2` (`k3`).
	*     - The slope at the end of the interval, using `k3` (`k4`).
	*   These increments are combined to produce the final value after the step.
	* - This method is popular for its accuracy and simplicity and is widely used in solving ODEs.
	* 
	* ## Memory Management
	* - The function uses `vector` to allocate temporary arrays (`dym`, `dyt`, `yt`) of size `n+1`.
	* - These arrays are freed at the end of the function using `free_vector`.
	* 
	* ## Notes
	* - The `y`, `dydx`, and `yout` arrays are indexed starting from 1 (not 0).
	* - The function signature uses non-standard C declaration style, which might be a holdover from older codebases or specific numerical libraries.
	*/

	int i;
	double xh,hh,h6,*dym,*dyt,*yt,*vector();
	void free_vector();

	dym=vector(1,n);
	dyt=vector(1,n);
	yt=vector(1,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=1;i<=n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_vector(yt,1,n);
	free_vector(dyt,1,n);
	free_vector(dym,1,n);
}
