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

typedef struct {
    //unlike with the system object above, this step_type object does not need
    //to match GSL's form explicitly, it just needs to define the method.
    //Since we are using RK-type methods exclusively, this is just a
    int dimension; //since we are passing a void pointer to do this, we need a way
    //to know how large it is in the end.
    void *butcher;
    //Make sure to put this at the end of the struct though,
    //in case we add more parts to it, nonspecific arrays must be the last element.
} nrpy_odiegm_step_type;