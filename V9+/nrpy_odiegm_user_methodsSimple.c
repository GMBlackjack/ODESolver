#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

//This file holds all the functions and definitions for the user to edit. 
//Note that it does not depend on any of the other files--so long as the formatting is maintained
//the operation of the code should be agnostic to what the user puts in here. 

//This structure here holds any constant parameters we may wish to report.
//Often this structure can be entirely empty if the equation is self-contained.
//But if we had a system that relied on an Equation of State, the parameters for that EOS would go here.

struct constantParameters { 
    int dimension; //number that says how many we have. 
    //add more as necessary.  
};

//Here's the prototypes     double rho;for the functions in this file, stated explicitly for the sake of clarity. 
void exceptionHandler (double x, double y[]); 
//handles any exceptions the user may wish to define.
int doWeTerminate (double x, double y[], struct constantParameters *params); //user-defined endpoint.
//Use if the code won't terminate itself from outside, or if there's a special condition. 
void constEval (double x, const double y[], struct constantParameters *params);
//Assign constants to the constantParameters struct based on values in y[]. 
int diffyQEval (double x, double y[], double dydx[], void *params);
//The definition for the system of equations itself goes here. 
int knownQEval (double x, double y[]);
//If an exact solution is known, it goes here, otherwise leave empty. 
void getInitialCondition (double y[]);
//Initial conditions for the system of differential equations. 
void assignConstants (double c[], struct constantParameters *params);
//Used to read values from constantParameters into an array so they can be reported in sequence. 

//Note that nrpy_odiegm_funcs.c does not depend on these definitions at all. The user is free
//to rename the functions if desired, though since diffyQEval and knownQEval are passed to 
//one of nrpy_odiegm's structs the actual function parameters for those two should not be adjusted.
//NOTE: the given nrpy_odiegm_main.c file will only work with the same names as listed here,
//only change names if creating a new custom main(). 

void exceptionHandler (double x, double y[])
{
    
}

int doWeTerminate (double x, double y[], struct constantParameters *params)
{
    if (x >= 2) {
        return 1;
    } else {
        return 0;
    }

}

void constEval (double x, const double y[], struct constantParameters *params)
{

}

int diffyQEval (double x, double y[], double dydx[], void *params)
{

        dydx[0] = y[0];
        //running the test on y'=y. Obviously gonna be e^x.

    return 1;
}


//This is the function to evaluate the known solution. Must be set manually.
int knownQEval (double x, double y[]) //This function is the other one passed using GSL's formulation. 
//Allows the specific_methods file to be completely agnostic to whatever the user is doing. 
{

    y[0] = exp(x);

    return 1;
    //report "success"
}

void getInitialCondition (double y[])
{
    y[0] = 1.0;
}

void assignConstants (double c[], struct constantParameters *params)
{

}