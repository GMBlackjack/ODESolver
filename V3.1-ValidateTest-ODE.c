#include <stdio.h>
#include <stdlib.h>
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
//This version is a "sideways" version that tests the order of the error term, it doesn't solve anything. 
//Code not guaranteed to make sense due to this.
//Very munimal user functionality, will need to be revamped. 
//Taylor Series Validation Method was removed since it was not appropriate for this context. 

//Heavily influenced by Numerical Mathematics and Computing 6E by Cheney and Kincaid.

//Outside the program, we substantiate the differential equation itself.
double diffyQEval (double x, double y)
{
    return y+1.0;
    //This is the differential equation itself. 
    //"return y" is the most basic result, for a basic exponential diffyQ.
    //feel free to change the return value to other functions. 
    //Note: not guaranteed to work for functions that are not well-behaved. 
}

//This is the function to evaluate the known solution. Must be set manually.
double knownQEval (double x)
{
    return exp(x)-1.0;
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
    double step = 0.0001; //the "step" value. 
    double bound = 0.0; //where the boundary/initial condition is.
    //Should work with non-integer boundaries.

    double bValue = 0.0; //the value at y(bound). 
    //by default we say y(0) = 1. 
    const int SIZE = 10000; //How many steps we are going to take. 
    bool validate = true; //set to true if you wish to run a validation test.
    //Attempts to find the order of the method used. 

    int method = 1; //sets the method. 1 is Euler's Method, 2 is Runge-Kutta Order 2. 4 is RK4
    //very basic user interface for method selection.
    printf("Which method to use? 1: Euler, 2: RK2, 4: RK4\n");
    scanf("%d", &method);
    printf("Using Method %d\n", method);


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

    double saveErr1=0, saveErr2=0; //variables for validation if requested. 
    
    //SECTION II: The Loop
    //prior to beginning the loop, start the timer. 

    double startN = time(NULL);
    //printf("Time: %.90f seconds\n", startN);
    //set start time to current time. Printing disabled.

    //also open the file we'll be writing data to. 
    FILE *fp;
    fp = fopen("oData.txt","w");

    //This loop fills out all the data. 
    //The switch is here to decide which method to use. 
    printf("INITIAL: Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound,yTruth1, y1, yError);
    fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound,yTruth1, y1, yError);
    switch(method) {
        case 1: {
            for (int i = 1; i < SIZE+1; i++){
        
                y2 = y1 + step*i*diffyQEval(bound,y1);
        
                yTruth2 = knownQEval(bound+step*i);
                yError = (yTruth2 - y2);

                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+i*step,yTruth2, y2, yError);
            }
            break;
            //Remember to use break or else switch keeps going. 
        }
        case 2: {
            double K1, K2; //The varabiles that store our Runge-Kutta results. 
            for (int i = 1; i < SIZE+1; i++){
                
                K1 = step*i*diffyQEval(bound,y1);
                K2 = step*i*diffyQEval(bound+step*i,y1 + K1);
                y2 = y1 + 0.5*(K1+K2);
        
                yTruth2 = knownQEval(bound+step*i);
                yError = (yTruth2 - y2);

                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+i*step,yTruth2, y2, yError);

            }
            break;
        }
        case 4: {
            double K1, K2, K3, K4; //The varabiles that store our Runge-Kutta results. 
            for (int i = 0; i < SIZE; i++){
                
                K1 = step*i*diffyQEval(bound,y1);
                K2 = step*i*diffyQEval(bound+0.5*step*i,y1 + K1*0.5);
                K3 = step*i*diffyQEval(bound+0.5*step*i,y1 + K2*0.5);
                K4 = step*i*diffyQEval(bound+step*i,y1 + K3);
                y2 = y1 + (1.0/6.0)*(K1+2.0*K2 + 2.0*K3 + K4);

                yTruth2 = knownQEval(bound+step*i);
                yError = (yTruth2 - y2);

                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+i*step,yTruth2, y2, yError);

            }
            break;
        }
        default: {
            printf("How did you get here? Did you use the wrong Method code?\n");
        }
    }

    //SECTION III: Analysis
    //Post-processing goes here.

    //DATA FORMATTER
    // basic reference: https://www.tutorialspoint.com/cprogramming/c_file_io.htm
    // used to be a file converter here, now there isn't, we just close the file. 
    fclose(fp);

    //TIMER
    double endN = time(NULL);
    //loop is complete, how long did it take?
    printf("Time Elapsed: %f seconds\n", endN-startN);
    //Only calculates to the closest second, for some reason.

    printf("ODE Solver \"Odie\" V1 Shutting Down...\n");
    return 0;
}

// - GM, master of dogs.