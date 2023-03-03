#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

struct constantParameters { 
    int dimension; //number that says how many we have. 
    double rho;
    //add more as necessary.  
};

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

int diffyQEval (double x, double y[], double dy[], void *params)
{

        dy[0] = y[0];
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

//Remember when adjusting these to adjust the necessary parameters in main() as well:
//step, bound, numberOfEquations, numberOfConstants, SIZE, and validate. 