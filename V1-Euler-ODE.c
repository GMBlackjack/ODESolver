#include <stdio.h>
#include <math.h>
#include <time.h> //we like to know how long things take. 
//having not used time.h before we refered to
//https://en.wikibooks.org/wiki/C_Programming/time.h

//Note: math.h requries the "-lm" arg be added at the END of tasks.json's arguments.
//https://askubuntu.com/questions/332884/how-to-compile-a-c-program-that-uses-math-h

//ODE Solver
//By G. M. Steward
//The main goal of this project is to solve Ordinary Differential Equation Systems
//in complete generality. However, we will take a while to get there.
//This first version, V1, has a simple goal: be able to solve first order
//ODE with a known result. 
//No user functionality yet, to adjust initial variables the code
//itself needs to be adjusted.

//Heavily influenced by Numerical Mathematics and Computing 6E by Cheney and Kincaid.

//Outside the program, we substantiate the differential equation itself.
double diffyQEval (double x, double y)
{
    return y;
    //This is the differential equation itself. 
    //"return y" is the most basic result, for a basic exponential diffyQ.
    //feel free to change the return value to other functions. 
    //Note: not guaranteed to work for functions that are not well-behaved. 
}

//This is the function to evaluate the known solution. Must be set manually.
double knownQEval (double x)
{
    return exp(x);
    //the known solution to the differential equaiton. 
    //used to measure relative errors. 
}

//Remember when adjusting these to adjust the boundary value bValue as well. 

int main()
{
    printf("Beginning ODE Solver \"Odie\" V1...\n");
    
    //SECTION I: Preliminaries
    //Before the program actually starts, variables need to be created
    //and set, as well as the function itself chosen. 
    //The diffyQ itself can be found declared in diffyQEval().
    double step = 0.001; //the "step" value. 
    int bound = 0; //where the boundary/initial condition is.
    //Functionality will have to be altered for non-integer boundaries.
    double bValue = 1; //the value at y(bound). 
    //by default we say y(0) = 1. 
    const int SIZE = 200000; //Size of the various arrays. 
    
    double y[SIZE];
    double yTruth[SIZE];
    double yError[SIZE];
    //These two arrays store funciton values. 
    //y is waht we solve for, yTruth is the "known results."
    //yError is the error of y when compared to yTruth
    //Make sure to set them to the same size. Default is 100 steps
    //Feel free to make large. 

    y[0] = bValue; //boundary condition, have to start somewhere.
    yTruth[0] = knownQEval(0);
    yError[0] = 0; //has to be zero as they must match at this point.
    
    printf("Position:\t%d\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\%\t\n",0,yTruth[0], y[0], yError[0]);
     
    //SECTION II: The Loop
    //prior to beginning the loop, start the timer. 

    double startN = time(NULL);

    printf("Time: %.90f seconds\n", startN);
    //set start time to current time. 

    //This loop will fill the arrays with information. 
    for (int i = 0; i < SIZE; i++){
        y[i+1] = y[i] + step*diffyQEval(i*step,y[i]);
        //This is Euler's method. 
        //Very basic, uses information on derivatives and functions
        //to calculate the function itself step by step.
        
        yTruth[i+1] = knownQEval(step*(i+1));
        yError[i+1] = (yTruth[i+1] - y[i+1])/yTruth[i+1];
        //Setting up the truth and error arrays at the same time.

        //After each step is calculated, print results. 
        printf("Position:\t%d\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\%\t\n",i+1,yTruth[i+1], y[i+1], yError[i+1]);
    }

    //SECTION III: Analysis
    //For now, there isn't much here, as the program prints its info as it find it. 
    //Post-processing will go here, which is why we've kept the array memory allocated, 
    // in case we want to use it here. 

    double endN = time(NULL);
    //loop is complete, how long did it take?
    printf("Time: %f seconds\n", startN);
    printf("Time: %f seconds\n", endN);
    printf("Time: %f seconds\n", endN-startN);
    //Only calculates to the closest second, for some reason.

    


    //Features:
    //Format Data for Excel for Visualization
    //Time of Calculation Analysis

    
    printf("ODE Solver \"Odie\" V1 Shutting Down...\n");
    return 0;
}

// - GM, master of version 1. 