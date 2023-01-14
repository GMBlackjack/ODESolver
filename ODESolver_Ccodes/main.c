#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdbool.h"
#include "time.h"

    //Outside the program, we substantiate the differential equation itself.
    double diffyQEval (double x, double y)
    {
        return exp(x) - exp(-x) - 1;
        //This is the differential equation itself. 
        //"return y+1.0" is the "default" selection for its simplicity yet usefulness for testing the algorithms.
        //feel free to change the return value to other functions. 
        //Note: not guaranteed to work for functions that are not well-behaved. 
    }

    //This is the function to evaluate the known solution. Must be set manually.
    double knownQEval (double x)
    {
        return exp(x) + exp(-x) - x;
        //the known solution to the differential equaiton. 
        //used to measure relative errors. 
        //"return exp(x)-1.0" is "default."
    }
    
/*
 * Fill out later
 */
int main() {

double butcher[3][3] = {{0.0,0,0},{1.0,1.0,0},{2.0,0.5,0.5}};
    printf("Beginning ODE Solver \"Odie\" V4...\n");
    
    //SECTION I: Preliminaries
    //Before the program actually starts, variables need to be created
    //and set, as well as the function itself chosen. 
    //The diffyQ itself can be found declared in diffyQEval().

    //How to get array size: https://stackoverflow.com/questions/37538/how-do-i-determine-the-size-of-my-array-in-c
    size_t methodSize = sizeof(butcher)/sizeof(butcher[0][0]);
    int dimension = sqrt((int)methodSize);
    printf("Method Order: %i. \nOrder of Error should be near Method Order + 1.\n",(int)butcher[dimension-1][0]);
    printf("If not, try a larger step size, roundoff error may be interfering.\n");

    double step = 0.001; //the "step" value.
    double bound = 0.0; //where the boundary/initial condition is.
    double bValue = 2.0; //the value at y(bound). By default we say y(0) = 0.
    const int SIZE = 1000; //How many steps we are going to take?
    bool validate = true; //set to true if you wish to run a validation test.
    //Attempts to find the order of the method used. 

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
    //If not, the boundary conditions were not stated properly. 

    double saveErr1=0, saveErr2=0; //variables for validation if requested. 
    
    //SECTION II: The Loop

    //prior to beginning the loop, start the timer. 
    double startN = time(NULL);
    //printf("Time: %.90f seconds\n", startN);
    //set start time to current time.  Uncomment to print.

    //also open the file we'll be writing data to. 
    FILE *fp;
    fp = fopen("oData.txt","w");

    //This loop fills out all the data.
    //It takes a provided butcher table and executes the method stored within. Any table should work. 
    printf("INITIAL: Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound,yTruth1, y1, yError);
    fprintf(fp,"Position:,\t%f,\tTruth:,\t%10.9e,\tCalculated:,\t%10.9e,\tError:,\t%10.9e,\t\n",bound,yTruth1, y1, yError);
    //Comma delimiters are printed to the file so it can be converted to .csv with ease. 
    for (int i = 0; i < SIZE; i++){
        
        double k[dimension]; //Create an array that can store all the steps needed for the method.
        for (int j = 1; j < dimension; j++) {
            //Due to the way the Butcher Table is formatted, start our index at 1 and stop at the end. 
            double xInsert = bound+i*step + butcher[j-1][0]*step;
            //x does not change much for different tables, just adjust the "step correction" term. 
            double yInsert = y1; //The formu of the y-insertion into our function changes, so it needs a loop.
            for (int n = 1; n < dimension; n++) {
                //Once again, start at index of 1 rather than 0. 
                yInsert = yInsert + butcher[j-1][n]*k[n];
            }
            k[j] = step*diffyQEval(xInsert,yInsert); //calculate the complete k-value. 
        }
        //Now that we have all the k-values set, we need to find the actual result in one final loop. 
        k[0] = y1;
        for (int j = 1; j < dimension; j++) {
            k[0] = k[0] + butcher[dimension-1][j]*k[j];
        }
        y2 = k[0];
        //possible efficiency option here: remove y2 and just use k[0]. That will make the code more obtuse, though.

        yTruth2 = knownQEval(bound+step*(i+1));
        yError = (yTruth2 - y2);

        //After each step is calculated, print results. 
        //printf("Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
        //uncomment if you want live updates. 
        fprintf(fp,"Position:,\t%f,\tTruth:,\t%10.9e,\tCalculated:,\t%10.9e,\tError:,\t%10.9e,\t\n",bound+(i+1)*step,yTruth2, y2, yError);
                
        //validation: grab the first nonzero error, calculate its order.
        if(validate==true && i == 0.0) {
            //Only activate on first step.
            saveErr1 = yError;

            //The following is an algorithm for determining the rate of error 
            //convergence. A bit rudimentary, could be condensed, but is also only
            //called once so not relaly a concern and it is easier to read this way. 
            for (int j = 1; j < dimension; j++) {
                double xInsert = bound+i*step*0.5 + butcher[j-1][0]*step*0.5;
                double yInsert = y1;
                for (int n = 1; n < dimension; n++) {
                    yInsert = yInsert + butcher[j-1][n]*k[n];
                }
                k[j] = step*0.5*diffyQEval(xInsert,yInsert);
            }
            k[0] = y1;
            for (int j = 1; j < dimension; j++) {
                k[0] = k[0] + butcher[dimension-1][j]*k[j];
            }
            double truthValidate = knownQEval(bound+step*0.5);
            saveErr2 = (truthValidate - k[0]);
            //Basically we just calculated the initial error for half step size. 
            //Now we can compare using the equation for order estimation:
            double order =  log2(saveErr1/saveErr2);
            printf("Order of Error: %f\n", order);
        }

        y1=y2;
        yTruth1=yTruth2;

    }

    //SECTION III: Analysis
    //Post-processing goes here.

    // basic reference: https://www.tutorialspoint.com/cprogramming/c_file_io.htm
    // used to be a file converter here, now there isn't, we just close the file. 
    fclose(fp);

    //TIMER
    double endN = time(NULL);
    //loop is complete, how long did it take?
    printf("Time Elapsed: %f seconds\n", endN-startN);
    //Only calculates to the closest second, for some reason.

    printf("ODE Solver \"Odie\" V4 Shutting Down...\n");
    return 0;
    }
