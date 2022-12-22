#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h> //we like to know how long things take. 
//having not used time.h before we refered to
//https://en.wikibooks.org/wiki/C_Programming/time.h

//Note: math.h requries the "-lm" arg be added at the END of tasks.json's arguments.
//https://askubuntu.com/questions/332884/how-to-compile-a-c-program-that-uses-math-h

//ODE Solver
//By G. M. Steward
//The main goal of this project is to solve Ordinary Differential Equation Systems
//in complete generality. However, we will take a while to get there.
//This second version, V2, seeks to add Runge-Kutta order 2
//Functionality for the solutions.
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
    //exp(x) is "default."
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
    double bound = 0; //where the boundary/initial condition is.
    //Functionality will have to be altered for non-integer boundaries.
    double bValue = 1; //the value at y(bound). 
    //by default we say y(0) = 1. 
    const int SIZE = 1000000; //How many steps we are going to take. 
    bool validate = true; //set to true if you wish to run a validation test

    int method = 2; //sets the method. 1 is Euler's Method, 2 is Runge-Kutta Order 2. 

    double y1;
    double y2;
    double yTruth1;
    double yTruth2;
    double yError;
    //These variables temporarily store the values calculated before they are 
    //printed to the output file and forgotten. 
    //y is what we solve for, yTruth is the "known results."
    //yError is the error of y when compared to yTruth
    //Should be able to handle any size, unlike arrays, which just hog memory. 

    y1 = bValue; //boundary condition, have to start somewhere.
    yTruth1 = knownQEval(bound);
    yError = 0; //has to be zero as they must match at this point.

    double saveErr1, saveErr2; //variables for validation if needed. 
    
    //SECTION II: The Loop
    //prior to beginning the loop, start the timer. 

    double startN = time(NULL);
    //printf("Time: %.90f seconds\n", startN);
    //set start time to current time. Printing disabled.

    //also open the file we'll be writing data to. 
    FILE *fp;
    fp = fopen("oData.txt","w");

    //This loop will fill the arrays with information. 
    //SIZE-1 required due to C's array indexing. 
    printf("INITIAL: Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound,yTruth1, y1, yError);
    fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound,yTruth1, y1, yError);
    switch(method) {
        case 1: {
            for (int i = 0; i < SIZE; i++){
        
                y2 = y1 + step*diffyQEval(bound+i*step,y1);
                //This is Euler's method. 
                //Very basic, uses information on derivatives and functions
                //to calculate the function itself step by step.
        
                yTruth2 = knownQEval(bound+step*(i+1));
                yError = (yTruth2 - y2)/yTruth2;
                //Setting up the truth and error arrays at the same time.

                //After each step is calculated, print results. 
                //printf("Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
                //uncomment if you want live updates. 
                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);

                if(validate==true && i == 0) {
                    saveErr1 = yError;
                    //validation requires two errors, we save it here. 
                    //we also need a half-step error. 
                    double temp = knownQEval(bound+step*0.5*(i+1));
                    saveErr2 = (temp - (y1 + step*0.5*diffyQEval(bound+i*step*0.5,y1)))/temp;

                }

                y1=y2;
                yTruth1=yTruth2;

            }
        }
        case 2: {
            double K1, K2; //The varabiles that store our Runge-Kutta results. 
            for (int i = 0; i < SIZE; i++){
                
                K1 = step*diffyQEval(bound+i*step,y1);
                K2 = step*diffyQEval(bound+i*step+step,y1 + K1);
                y2 = y1 + step*diffyQEval(bound+i*step,y1);
                //This is the Runge-Kutta 2 method.  
                //Should have second-order error. 
        
                yTruth2 = knownQEval(bound+step*(i+1));
                yError = (yTruth2 - y2)/yTruth2;
                //Setting up the truth and error arrays at the same time.

                //After each step is calculated, print results. 
                //printf("Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
                //uncomment if you want live updates. 
                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);

                if(validate==true && i == 0) {
                    saveErr1 = yError;
                    //validation requires two errors, we save it here. 
                    //we also need a half-step error. 
                    double temp = knownQEval(bound+step*0.5*(i+1));
                    saveErr2 = (temp - (y1 + step*0.5*diffyQEval(bound+i*step*0.5,y1)))/temp;

                }

                y1=y2;
                yTruth1=yTruth2;

            }
        }
    }

    //SECTION III: Analysis
    //For now, there isn't much here, as the program prints its info as it finds it. 
    //Post-processing will go here.

    //DATA FORMATTER
    // basic reference: https://www.tutorialspoint.com/cprogramming/c_file_io.htm
    // used to be a file converter here, now there isn't, we just close the file. 
    fclose(fp);

    //TIMER
    double endN = time(NULL);
    //loop is complete, how long did it take?
    printf("Time Elapsed: %f seconds\n", endN-startN);
    //Only calculates to the closest second, for some reason.

    //validation test, only executes if requested
    //also very quick and dirty, doesn't work for all functions, needs revamp. 
    if (validate == true){
        printf("Validating...\n");
        
        double yError2;

        y1 = bValue; 
        yTruth1 = knownQEval(bound);
        yError = 0;
        //Resetting everything for easy access 

        //we only need two things, so we can reuse y1 and y2s values.
        y2 = y1 + step*diffyQEval(bound,y1);
        yTruth2 = knownQEval(bound+step);
        yError = (yTruth2 - y2)/yTruth2;
        
        y1 = y1 + step*0.5*diffyQEval(bound,y1);
        yTruth1 = knownQEval(bound+step*0.5);
        yError2 = ((yTruth1 - y1)/yTruth1);

        //we now have two first-step approximations, one according to the bound and another according to
        //half the step.

        double order =  log2(yError/yError2);
        printf("Order: %f\n", order);
        //this notably only works for certain functions. For instance, y' = x/(y+1), y(0)=0 always produces the exact
        //same error (100%) on the first step since all the initial evaluations are zero every time. 
    }

    printf("ODE Solver \"Odie\" V1 Shutting Down...\n");
    return 0;
}

// - GM, master of version 1. 