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
//This second version, V3, is mainly concerned with validations.
//As well as comment cleanup. 
//Very munimal user functionality, will need to be revamped. 

//Heavily influenced by Numerical Mathematics and Computing 6E by Cheney and Kincaid.

//Outside the program, we substantiate the differential equation itself.
double diffyQEval (double x, double y)
{
    return x/(y+1.0);
    //This is the differential equation itself. 
    //"return y" is the most basic result, for a basic exponential diffyQ.
    //feel free to change the return value to other functions. 
    //Note: not guaranteed to work for functions that are not well-behaved. 
}

//This is the function to evaluate the known solution. Must be set manually.
double knownQEval (double x)
{
    return sqrt(x*x+1.0)-1.0;
    //the known solution to the differential equaiton. 
    //used to measure relative errors. 
    //exp(x) is "default."
}

//This funciton evaluates the Taylor series expansion for e^x
//Used for validation. Not helpful if you adjusted the functions above.
//Completely optoinal though. 
double taylorQEval (double x, double order) {
    double result = 0;
    if (order > 0) {
        result = result + 0;
    }
    if (order > 1) {
        result = result + (x*x)/2.0;
    }
    if (order > 2) {
        result = result + 0;
    }
    if (order > 3) {
        result = result - (x*x*x*x)/8.0;
    }
    if (order > 4) {
        result = result + 0;
    }
    return result;
}

//Remember when adjusting these to adjust the boundary value bValue as well. 

int main()
{
    printf("Beginning ODE Solver \"Odie\" V1...\n");
    
    //SECTION I: Preliminaries
    //Before the program actually starts, variables need to be created
    //and set, as well as the function itself chosen. 
    //The diffyQ itself can be found declared in diffyQEval().
    double step = 0.01; //the "step" value. 
    double bound = 0; //where the boundary/initial condition is.
    //Functionality will have to be altered for non-integer boundaries.
    double bValue = 0; //the value at y(bound). 
    //by default we say y(0) = 1. 
    const int SIZE = 100; //How many steps we are going to take. 
    bool validate = true; //set to true if you wish to run a validation test

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

    double saveErr1=0, saveErr2=0, saveIndex = 0; //variables for validation if needed. 
    //validation is currently under construction. 
    
    //SECTION II: The Loop
    //prior to beginning the loop, start the timer. 

    double startN = time(NULL);
    //printf("Time: %.90f seconds\n", startN);
    //set start time to current time. Printing disabled.

    //also open the file we'll be writing data to. 
    FILE *fp;
    fp = fopen("oData.txt","w");

    //This loop fills out all the data. The containing switch is there to decide which method to use. 
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
                yError = (yTruth2 - y2);

                //After each step is calculated, print results. 
                printf("Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
                //uncomment if you want live updates. 
                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
                
                //validation: grab the first nonzero error, save its location.
                if(validate==true && sqrt(saveErr1*saveErr1) <= 0.000000000000001) {
                    //tiny but nonzero number allows us to avoid roundoff error. 
                    saveErr1 = yError;
                    saveIndex = bound+step*(i+1);

                    double yValidate = y1 + step*0.5*diffyQEval(bound+i*step*0.5,y1);
                    double truthValidate = knownQEval(bound+step*0.5);
                    saveErr2 = (truthValidate - yValidate);
                    double order =  log2(saveErr1/saveErr2);
                    printf("Order of Error: %f\n", order);
                }

                y1=y2;
                yTruth1=yTruth2;

            }
            break;
            //Remember to use break or else switch keeps going. 
        }
        case 2: {
            double K1, K2; //The varabiles that store our Runge-Kutta results. 
            for (int i = 0; i < SIZE; i++){
                
                K1 = step*diffyQEval(bound+i*step,y1);
                K2 = step*diffyQEval(bound+i*step+step,y1 + K1);
                y2 = y1 + 0.5*(K1+K2);
                //This is the Runge-Kutta 2 method.  
                //Should have second-order error. 
        
                yTruth2 = knownQEval(bound+step*(i+1));
                yError = (yTruth2 - y2);

                //After each step is calculated, print results. 
                printf("Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
                //uncomment if you want live updates. 
                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);

                //validation: grab the first nonzero error, save its location.
                //validation: grab the first nonzero error, save its location.
                if(validate==true && sqrt(saveErr1*saveErr1) <= 0.000000000000001) {
                    //tiny but nonzero number allows us to avoid roundoff error. 
                    saveErr1 = yError;
                    saveIndex = bound+step*(i+1);

                    K1 = step*0.5*diffyQEval(bound+i*step,y1);
                    K2 = step*0.5*diffyQEval(bound+i*step*0.5+step*0.5,y1 + K1);
                    double yValidate = y1 + 0.5*(K1+K2);
                    double truthValidate = knownQEval(bound+step*0.5);
                    saveErr2 = (truthValidate - yValidate);
                    double order =  log2(saveErr1/saveErr2);
                    printf("Order of Error: %f\n", order);
                }

                y1=y2;
                yTruth1=yTruth2;

            }
            break;
        }
        case 4: {
            double K1, K2, K3, K4; //The varabiles that store our Runge-Kutta results. 
            for (int i = 0; i < SIZE; i++){
                
                K1 = step*diffyQEval(bound+i*step,y1);
                K2 = step*diffyQEval(bound+i*step+0.5*step,y1 + K1*0.5);
                K3 = step*diffyQEval(bound+i*step+0.5*step,y1 + K2*0.5);
                K4 = step*diffyQEval(bound+i*step+step,y1 + K3);
                y2 = y1 + (1.0/6.0)*(K1+2.0*K2 + 2.0*K3 + K4);
                //This is the Runge-Kutta 4 method.  
                //Should have fourth-order error. 
                //Notably gives answers reasonably accurate up to rounding error! 

                yTruth2 = knownQEval(bound+step*(i+1));
                yError = (yTruth2 - y2);

                //After each step is calculated, print results. 
                printf("Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
                //uncomment if you want live updates. 
                fprintf(fp,"Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);

                //validation: grab the first nonzero error, save its location.
                if(validate==true && sqrt(saveErr1*saveErr1) <= 0.000000000000001) {
                    //tiny but nonzero number allows us to avoid roundoff error. 
                    saveErr1 = yError;
                    saveIndex = bound+step*(i+1);

                    K1 = step*0.5*diffyQEval(bound+i*step*0.5,y1);
                    K2 = step*0.5*diffyQEval(bound+i*step*0.5+0.25*step,y1 + K1*0.5);
                    K3 = step*0.5*diffyQEval(bound+i*step*0.5+0.25*step,y1 + K2*0.5);
                    K4 = step*0.5*diffyQEval(bound+i*step*0.5+0.5*step,y1 + K3);
                    double yValidate = y1 + (1.0/6.0)*(K1+2.0*K2 + 2.0*K3 + K4);
                    double truthValidate = knownQEval(bound+step*0.5);
                    saveErr2 = (truthValidate - yValidate);
                    double order =  log2(saveErr1/saveErr2);
                    printf("Order of Error: %f\n", order);
                }

                y1=y2;
                yTruth1=yTruth2;

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

    //validation test, only executes if requested
    if (validate == true){
        printf("Validating...\n");
        
        double taylorVal;

        taylorVal = taylorQEval(saveIndex, method);
        yTruth1 = knownQEval(saveIndex);
        yError = (yTruth1 - taylorVal);

        printf("Index: %f Numerical: %10.9e - Taylor: %10.9e = %10.9e\n",saveIndex, saveErr1, yError, saveErr1-yError);
        printf("If result is near zero (1e15 or smaller is sufficient) then the order is correct.\n");
    }

    printf("ODE Solver \"Odie\" V1 Shutting Down...\n");
    return 0;
}

// - GM, master of version 1. 