#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

//ETK interface. 
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// This file holds all the functions and definitions for the user to edit. 
// Note that it does not depend on any of the other files--so long as the formatting is maintained
// the operation of the code should be agnostic to what the user puts in here. 

// This struct here holds any constant parameters we may wish to report.
// Often this struct can be entirely empty if the system of equations is self-contained.
// But if we had a system that relied on an Equation of State, 
// the parameters for that EOS would go here. 
struct constant_parameters { 
    int dimension; // number that says how many constants we have. 
    double rho;
    // double parameter;
    // add more as necessary. Label as desired. 
};

// Here are the prototypes for the functions in this file, stated explicitly for the sake of clarity. 
void exception_handler (double x, double y[]); 
// Handles any exceptions the user may wish to define.
int do_we_terminate (double x, double y[], struct constant_parameters *params); 
// User-defined endpoint.
// Generally used if the code won't terminate itself from outside, or if there's a variable condition. 
void const_eval (double x, const double y[], struct constant_parameters *params);
// Assign constants to the constant_parameters struct based on values in y[]. 
int diffy_Q_eval (double x, double y[], double dydx[], void *params);
// The definition for the system of equations itself goes here. 
int known_Q_eval (double x, double y[]);
// If an exact solution is known, it goes here, otherwise leave empty. 
void get_initial_condition (double y[]);
// Initial conditions for the system of differential equations. 
void assign_constants (double c[], struct constant_parameters *params);
// Used to read values from constant_parameters into an array so they can be reported in sequence. 

// Note that nrpy_odiegm_funcs.c does not depend on these definitions at all. The user is free
// to rename the functions if desired, though since diffy_Q_eval and known_Q_eval are passed to 
// one of nrpy_odiegm's structs the actual function parameters for those two should not be adjusted.
// NOTE: the given nrpy_odiegm_main.c file will only work with the same names as listed here,
// only change names if creating a new custom main function. 

void exception_handler (double x, double y[])
{
    // This funciton might be empty. It's only used if the user wants to hard code some limitations 
    // on some varaibles.
    // Good for avoding some divide by zero errors, or going negative in a square root. 
    if (y[0] < 0) {
        y[0] = 0;
    }
    // In this case, the TOV Equations, we need to make sure the pressure doesn't go negative.
    // Physically, it cannot, but approximation methods can cross the P=0 line
    // We just need a hard wall to prevent that. 
}

int do_we_terminate (double x, double y[], struct constant_parameters *params)
{
    DECLARE_CCTK_PARAMETERS
    // This funciton might be empty. It's only used if the user wants to have 
    // a special termination condition.
    // Today we do. We terminate once the pressure hits zero, or goes below it. 
    // Notably we also consider ridiculously small pressures to be "zero" since we might be asymptotic. 
    if (x > TOVOdieGM_outer_radius_limit) {
        return 1;
    } else {
        return 0;
    }
    // return 1; for termination.
}

void const_eval (double x, const double y[], struct constant_parameters *params)
{
    DECLARE_CCTK_PARAMETERS
    // Sometimes we want to evaluate constants in the equation that change, 
    // but do not have derivative forms.
    // Today, we do that for the total energy density. 
    // params->rho = sqrt(y[0]) + y[0];
    params->rho = pow(y[0] / TOVOdieGM_K , 1.0 / TOVOdieGM_Gamma) + y[0] / (TOVOdieGM_Gamma - 1.0);
    // The total energy density only depends on pressure. 
}

int diffy_Q_eval (double x, double y[], double dydx[], void *params)
{
    // GSL-adapted evaluation function. 
    // It is possible to do this with one array, but GSL expects two. 

    // Always check for exceptions first, then perform evaluations. 
    exception_handler(x,y);
    const_eval(x,y,params);

    // Dereference the struct
    double rho = (*(struct constant_parameters*)params).rho;
    // double parameter = (*(struct constant_parameters*)params).parameter;
    // WHY oh WHY GSL do you demand we use a VOID POINTER to the struct...?
    // https://stackoverflow.com/questions/51052314/access-variables-in-struct-from-void-pointer
    // Make sure to dereference every parameter within the struct so it can be used below. 

    // This if statement is an example of a special condition, 
    // in this case at x=0 we have a divide by zero problem. 
    // Fortunately, we manually know what the derivatives should be.
    // Alternatively, we could define piecewise equations this way. 
    if(x == 0) {
        dydx[0] = 0; 
        dydx[1] = 0;
        dydx[2] = 0;
        dydx[3] = 1;
    }
    else {
        dydx[0] = -((rho+y[0])*( (2.0*y[2])/(x) + 8.0*3.1415926535897931160*x*x*y[0] ))/(x*2.0*(1.0 - (2.0*y[2])/(x)));
        dydx[1] =  ((2.0*y[2])/(x) + 8.0*3.1415926535897931160*x*x*y[0])/(x*(1.0 - (2.0*y[2])/(x)));
        dydx[2] = 4*3.1415926535897931160*x*x*rho;
        dydx[3] = (y[3])/(x*sqrt(1.0-(2.0*y[2])/x));
        // Visual Studio likes to complain that M_PI is not defined, even though it is. 
        // So we used 3.1415926535897931160. which is just M_PI printed out to extra digits.
        // There was no observed change in the final product. 
    }
    // This funciton is not guaranteed to work in all cases. For instance, we have manually 
    // made an exception for x=0, since evaluating at 0 produces infinities and NaNs. 
    // Be sure to declare any exceptions before running, both here and in exception_handler, 
    // depending on the kind of exception desired.  

    return 0;
    // GSL_SUCCESS is 0. We do not support fancy error codes like GSL. 
}

// This is the function to evaluate the known solution. Must be set manually.
int known_Q_eval (double x, double y[]) // This function is another one passed using GSL's formulation. 
// Allows the nrpy_odiegm_user_methods.c file to be completely agnostic to whatever the user is doing. 
{
    // y[0] = ...
    // y[1] = ...
    // This function is only used if there are known solutions. 
    // Notably this is not the case for the TOV equations. 
    // If you do put anything here, make SURE it has the same order as the differential equations. 
    // In the case of TOV, that would be Pressure, nu, mass, and r-bar, in that order. 

    return 1;
    // report "success," what would have been GSL_SUCCESS in the GSL formulation. 
}

void get_initial_condition (double y[])
{
    DECLARE_CCTK_PARAMETERS
    // be sure to have these MATCH the equations in diffy_Q_eval
    y[0] = TOVOdieGM_K*pow(TOVOdieGM_central_baryon_density,TOVOdieGM_Gamma); // Pressure, can be calcualated from central baryon density. 
    y[1] = 0.0; // nu
    y[2] = 0.0; // mass
    y[3] = 0.0; // r-bar
}

void assign_constants (double c[], struct constant_parameters *params)
{
    // Reading parameters from the constant_parameters struct is rather difficult, since it exists
    // in the higher order "objects" as a void pointer. So the user should declare what constants
    // are what for ease of use, usually for printing in an algorithmic way.
    c[0] = params->rho; // Total energy density. 
    // Add more as required. 
}
