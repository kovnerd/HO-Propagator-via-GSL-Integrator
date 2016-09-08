#include <stdlib.h>
#include <math.h>
#define TIME 4.
#define N 8
#define LATTICE_SPACING 0.5 // T/N
#define MASS 1.

double lattice_action(double, size_t, double *);
double integrand(double *, size_t, void*);
double potential(double);

