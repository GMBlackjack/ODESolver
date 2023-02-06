#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

typedef struct {
    int (*function) (double x, const double y[], double dydx[], void *params);
    //The function passed to this struct contains the definitions of the differnetial equations. 
    int (*jacobian) (double t, const double y[], double *dfdy, double dfdt[],
                   void *params); //The jacobian is a holdover from GSL, it will not be used in this program,
                   //Always pass NULL to this part. 
                   //If we were not making a drop-in replacement, we would just have a three-paraemeter struct.
    size_t dimension; //For storing how big our system of equations is. 
    //Just pass it an int, usually. 
    void *params; //For storing extra constants needed to evaluate the functions. 
} nrpy_odiegm_system;

