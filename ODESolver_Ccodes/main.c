#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdbool.h"
#include "time.h"

    double diffyQEval (double x, double y[], int i)
    {
        switch(i){
            case 0: {
                return y[1]; 
                //we use 0 to represent the answer we eventually want to get
                //in this case, the solution to y'' = y+x, stated here as y'=z
                break;
            }
            case 1: {
                return y[0] + x;
                //1 represents the first derivative. y''=y+x, but when we split it up we got y'=z. This equation returns the value
                //for the z derivative, z'=y+x
                break;
            }
            // case...
            //add more cases for higher order systems. 
            default: {
                //Any number could be a result, theoretically, so there isn't an easy error value to put here.
            }
        }
        //This is the differential equation system itself. 
        //By default we have a very simple y'' = y+x situation here, split up into
        // y[0]' = y[1]
        // y[1]' = y[0]+x
        //Naturally other equaitons can be put in, but be sure to change the numberOfEquations value!
        //feel free to change the return values to other functions. 
        //Note: not guaranteed to work for functions that are not well-behaved. 
        //We have this set up so y[0] would be the final "answer" for a split-up higher-order ODE. 
    }

    //This is the function to evaluate the known solution. Must be set manually.
    double knownQEval (double x)
    {
        return exp(x) + exp(-x) - x;
        //the known solution to the differential equaiton, specifically what we call y[0]
        //used to measure relative errors. 
        //"exp(x) + exp(-x) - x" is "default." as it is the answer to the differential equation we chose. 
        //Do note that this would change with different boundary conditions. 
    }
    
/*
 * Fill out later
 */
int main() {

double butcher[5][5] = {{0.0,0,0,0,0},{0.5,0.5,0,0,0},{0.5,0.0,0.5,0,0},{1.0,0.0,0.0,1.0,0},{4.0,0.16666666666666666,0.3333333333333333,0.3333333333333333,0.16666666666666666}};
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

    double step = 0.1; //the "step" value.
    double bound = 0.0; //where the boundary/initial condition is.
    int numberOfEquations = 2; //How many equations are in our system?
    //Be very careful setting these boundary conditions, they need to match the number of equations. 
    double bValue[numberOfEquations]; 
    bValue[0] = 2.0;
    bValue[1] = -1.0;
    //the value at y(bound). By default we say y(0) = 0.
    const int SIZE = 10; //How many steps we are going to take?
    bool validate = true; //set to true if you wish to run a validation test.
    //Attempts to find the order of the method used. 

    double y1[numberOfEquations];
    double y2[numberOfEquations];
    double yTruth1;
    double yTruth2;
    double yError;
    //These variables temporarily store the values calculated before they are 
    //printed to the output file and forgotten. 
    //y is what we solve for, yTruth is the "known results."
    //yError is the error of y when compared to yTruth
    //Should be able to handle any size, unlike arrays, which just hog memory. 
    //Errors and Truth are not arrays as they are only concerned with the final value. 

    for (int n = 0; n < numberOfEquations; n++) {
        y1[n] = bValue[n];
        //Assign the initial values to the boundary conditions. 
    } 
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
    printf("INITIAL: Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound,yTruth1, y1[0], yError);
    fprintf(fp,"Position:,\t%f,\tTruth:,\t%10.9e,\tCalculated:,\t%10.9e,\tError:,\t%10.9e,\t\n",bound,yTruth1, y1[0], yError);
    //Comma delimiters are printed to the file so it can be converted to .csv with ease. 
    
    for (int i = 0; i < SIZE; i++){ 
        //printf("%i",m);
        double K[dimension][numberOfEquations];
        //Since we'll be calling this while it's empty, even though there should be no errors due
        //to the way it's set up, let's go ahead and fill it with zeroes.
        for (int j = 0; j<dimension; j++) {
            for (int n = 0; n<numberOfEquations; n++) {
                K[j][n]=0.0;
            }
        }
        //each diffyQ has its own set of K-values, one for each equation. 

        double yInsert[numberOfEquations];
        //We also need an array for the inserted y-values for each equation. 
        //Most applications actually have the different yInsert values be independent, so 
        //if we knew the form of the equation we could simplify the code.
        //However, we need to make sure to always fill everything in case we have a system
        //of the form y'=f(u,y) u'=g(u,y)

        for (int j = 1; j < dimension; j++) {
            //Due to the way the Butcher Table is formatted, start our index at 1 and stop at the end. 
            double xInsert = bound+i*step + butcher[j-1][0]*step;
            //x does not change much for different tables, just adjust the "step correction" term.
            //Is the same for every equation too.

            for (int n = 0; n < numberOfEquations; n++) {
                yInsert[n] = y1[n];
            } 

            for (int n = 1; n < dimension; n++) {
                //Once again, start at index of 1 rather than 0.
                for (int q = 0; q < numberOfEquations; q++) {
                    yInsert[q] = yInsert[q] + butcher[j-1][n]*K[n][q];
                }
                //Each individual y portion is dependent on one of the K values.
                //K values are initially set to zero even though technically whenever 
                //we would use an undeclared K-value the butcher table would have zero.
                //You know, just in case something goes wrong. 
            }

            for (int n = 0; n < numberOfEquations; n++) {
                K[j][n] = step*diffyQEval(xInsert,yInsert,n);
                //Actually calculate the K-values.
            } 

        }
        //Now that we have all the K-values set, we need to find the actual result in one final loop.
        //The sum for the first set... 
        for (int n = 0; n< numberOfEquations; n++) {
            K[0][n] = y1[n];
            for (int j = 1; j < dimension; j++) {
                K[0][n] = K[0][n] + butcher[dimension-1][j]*K[j][n];
            }
            y2[n] = K[0][n];
        }
        
        yTruth2 = knownQEval(bound+step*(i+1));
        yError = (yTruth2 - y2[0]);

        //After each step is calculated, print results. 
        //printf("Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound+(i+1)*step,yTruth2, y2, yError);
        //uncomment if you want live updates. 
        //if (m==0) {
            fprintf(fp,"Position:,\t%f,\tTruth:,\t%10.9e,\tCalculated:,\t%10.9e,\tError:,\t%10.9e,\t\n",bound+(i+1)*step,yTruth2, y2[0], yError);
        //}
                
        //validation: grab the first nonzero error, calculate its order.
        //Currently broken. 
        if(validate==true && i == 0.0) { //currently set to ignore this. REMEMBER TO TURN BACK ON LATER!
            //Only activate on first step. 
            saveErr1 = yError;

            //The following is an algorithm for determining the rate of error 
            //convergence. A bit rudimentary, could be condensed, but is also only
            //called once so not relaly a concern and it is easier to read this way. 
            double step2 = step*0.5;
            //It is easier to just use another variable than multiply the step by 0.5 every time. 
                for (int j = 0; j<dimension; j++) {
                    for (int n = 0; n<numberOfEquations; n++) {
                        K[j][n]=0.0;
                    }
                }
                for (int j = 1; j < dimension; j++) {
                    double xInsert = bound+i*step2 + butcher[j-1][0]*step2;
                    for (int n = 0; n < numberOfEquations; n++) {
                        yInsert[n] = y1[n];
                    } 
                    for (int n = 1; n < dimension; n++) {
                        for (int q = 0; q < numberOfEquations; q++) {
                            yInsert[q] = yInsert[q] + butcher[j-1][n]*K[n][q];
                        }
                    }
                    for (int n = 0; n < numberOfEquations; n++) {
                        K[j][n] = step2*diffyQEval(xInsert,yInsert,n);
                    } 
                }
                for (int n = 0; n< numberOfEquations; n++) {
                    K[0][n] = y1[n];
                    for (int j = 1; j < dimension; j++) {
                        K[0][n] = K[0][n] + butcher[dimension-1][j]*K[j][n];
                    }
                }
            double truthValidate = knownQEval(bound+step2);
            saveErr2 = (truthValidate - K[0][0]);
            //Basically we just calculated the initial error for half step size. 
            //Now we can compare using the equation for order estimation:
            double order =  log2(saveErr1/saveErr2);
            printf("Order of Error: %f\n", order);
        }
        
        for (int n = 0; n < numberOfEquations; n++) {
            y1[n]=y2[n];
        } 
        //make sure to assign all variables to the next step. 
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

// - GM, master of dogs.
    }
