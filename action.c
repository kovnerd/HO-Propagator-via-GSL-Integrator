#include "action.h"

/*
 * NOTE: In this implementation, the parameter "xs" to the lattice_action() and integrand() functions represent 
 * the possible path integral paths:
 * 		X_1, X_2, ... X_N-1
 * 
 * this parameter does NOT include information about X_0 or X_N, as in this implementations, the boundary conditions are:
 * 
 * 		X_0 = X_N
 * 
 */ 
 
double lattice_action(double x_0, size_t d, double *xs)//SOMETHING ISN'T RIGHT HERE
{
	double xNow2 = 0.0, xNext2 = 0.0, xNowNext = 0.0;//quadratic term of the lattice action
	double v_term = 0.;//potential-dependent term in the action
	int dim = (int)d;
	for(int i = 0; i < N; i++)//sum from 0 to N - 1
	{
		//implementation of boundary conditions
		double xNext = (i == dim) ? x_0 : xs[i];
		double xNow = (i == 0) ? x_0 : xs[i - 1];
		
		//calculation of lattic action terms
		
		xNow2 += xNow * xNow;
		xNext2 += xNext * xNext;
		xNowNext += 2 * xNext * xNow;
		v_term += potential(xNow);
	}
	
	return ((0.5*MASS)/LATTICE_SPACING) * (xNow2 + xNext2 - xNowNext) + (LATTICE_SPACING * v_term);
}

double integrand(double *xs, size_t dim, void *params)
{
	double x_0 = *((double *) params);
	double S_lat = lattice_action(x_0, dim, xs);
	
	return exp(-S_lat);
}

double potential(double x)
{
	return 0.5*MASS*x*x*x*x;
}
