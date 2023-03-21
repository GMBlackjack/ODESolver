#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

//This file holds all the functions and definitions for the user to edit. 

//This structure here holds any constant parameters we may wish to report.
//Often this structure can be entirely empty if the equation is self-contained.
//But if we had a system that relied on an Equation of State, the parameters for that EOS would go here. 
struct constantParameters { 
    int dimension; //number that says how many we have. 
    double rho;
    //double parameter;
    //add more as necessary. Label as desired. 
};

void exceptionHandler (double x, double y[])
{
    //This funciton might be empty. It's only used if the user wants to hard code some limitations 
    //On some varaibles.
    //Good for avoding some divide by zero errors, or going negative in a square root. 
    if (y[0] < 0) {
        y[0] = 0;
    }
    //In this case, the TOV Equations, we need to make sure the pressure doesn't go negative.
    //Physically, it cannot, but approximation methods can cross the P=0 line
    //We just need a hard wall to prevent that. 
}

int doWeTerminate (double x, double y[], struct constantParameters *params)
{
    //This funciton might be empty. It's only used if the user wants to have a special termination condition.
    //Today we do. We terminate once the pressure hits zero, or goes below it. 
    //Notably we also consider ridiculously small pressures to be "zero" since we might be asymptotic. 
    if (y[0] < 1e-16) {
        return 1;
    } else {
        return 0;
    }
    //return 1 for termination.
}

void constEval (double x, const double y[], struct constantParameters *params)
{
    //Sometimes we want to evaluate constants in the equation that change, 
    //but do not have derivative forms.
    //Today, we do that for the total energy density. 
    params->rho = sqrt(y[0]) + y[0];
    //The total energy density only depends on pressure. 
}

int diffyQEval (double x, double y[], double dydx[], void *params)
{
    //GSL-adapted evaluation function. 
    //It is possible to do this with one array, but GSL expects two. 

    //Always check for exceptions first, then perform evaluations. 
    exceptionHandler(x,y);
    constEval(x,y,params);

    //dereference the struct
    double rho = (*(struct constantParameters*)params).rho;
    //double parameter = (*(struct constantParameters*)params).parameter;
    //WHY oh WHY GSL do you demand we use a VOID POINTER to the struct...
    //https://stackoverflow.com/questions/51052314/access-variables-in-struct-from-void-pointer
    //make sure to dereference every parameter within the struct so it can be used below. 

    //This if statement is an example of a special condition, 
    //in this case at x=0 we have a divide by zero problem. 
    //Fortunately, we manually know what the derivatives should be.
    //Alternatively, we could define piecewise equations this way. 
    if(x == 0) {
        dydx[0] = 0; 
        dydx[1] = 0;
        dydx[2] = 0;
        dydx[3] = 1;
    }
    else {
        dydx[0] = -((rho+y[0])*( (2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0] ))/(x*2.0*(1.0 - (2.0*y[2])/(x)));
        dydx[1] =  ((2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0])/(x*(1.0 - (2.0*y[2])/(x)));
        dydx[2] = 4*M_PI*x*x*rho;
        dydx[3] = (y[3])/(x*sqrt(1.0-(2.0*y[2])/x));
        //Visual Studio likes to complain that M_PI is not defined, even though it is. 
    }
    //This funciton is not guaranteed to work in all cases. For instance, we have manually 
    //made an exception for x=0, since evaluating at 0 produces infinities and NaNs. 
    //Be sure to declare any exceptions before running, both here and in exceptionHandler(), depending 
    //on the kind of exception desired.  

    return 0;
    //GSL_SUCCESS is 0. We do not support fancy error codes like GSL. 
}

//This is the function to evaluate the known solution. Must be set manually.
int knownQEval (double x, double y[]) //This function is another one passed using GSL's formulation. 
//Allows the specific_methods file to be completely agnostic to whatever the user is doing. 
{
    //y[0] = ...
    //y[1] = ...
    //This function is only used if there are known solutions. 
    //Notably this is not the case for the TOV equations. 
    //If you do put anything here, make SURE it has the same order as the differential equations. 
    //In the case of TOV, that would be Pressure, nu, mass, and r-bar, in that order. 

    return 1;
    //report "success," what would have been GSL_SUCCESS in the original formulation. 
}

void getInitialCondition (double y[])
{
    //be sure to have these MATCH the equations in diffyQEval
    y[0] = 0.016714611225000002; //Pressure, can be calcualated from central baryon density. 
    y[1] = 0.0; //nu
    y[2] = 0.0; //mass
    y[3] = 0.0; //r-bar
}

void assignConstants (double c[], struct constantParameters *params)
{
    //Reading parameters from the constantParameters struct is rather difficult, since it exists
    //in the higher order "objects" as a void pointer. So the user needs to declare what constants
    //are what. 
    c[0] = params->rho; //total energy density. 
    //add more as required. 
}

//Remember when adjusting these to adjust the necessary parameters in main() as well:
//step, bound, numberOfEquations, numberOfConstants, SIZE, and validate. 