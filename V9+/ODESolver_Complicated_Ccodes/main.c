#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdbool.h"
#include "time.h"

void diffyQEval (double x, double y[], double c[])
{
    //Efficient(ish) Assignment Method
    //This computes the values of the *derivatives* of the functions, so what we are passed
    //in the array y[] is *not* the same thing as what we are calculating, be careful not to re-use
    //outside of the function, because we are taking in the function itself and replacing it with its derivative. 
    //This is somewhat counter-intuitive, but it is efficient. 

    //This if statement is an example of a special condition, in this case at x=0 we have a divide by zero problem. 
    //In this case we manually know what the derivatives should be.
    //Alternatively, we could define piecewise equations this way. 
    if(x == 0) {
        y[0] = 0; 
        y[1] = 0;
        y[2] = 0;
        y[3] = 1;
    }
    else {
        double y0 = y[0]; //Pressure
        double y1 = y[1]; //nu. Uncertain what this physically is, likely just a bookkeeping situation. 
        double y2 = y[2]; //Total mass. 
        double y3 = y[3]; //r-bar, the isotropic radius. NOT NORMALIZED. Must be normalized in post-processing. 

        //Note that we declare these buffer variables since we will need to use them multiple times, but we are also
        //overwriting them in the original array, so we need buffers. 
        y[0] = -((c[0]+y0)*( (2.0*y2)/(x) + 8.0*M_PI*x*x*y0 ))/(x*2.0*(1.0 - (2.0*y2)/(x)));
        y[1] =  ((2.0*y2)/(x) + 8.0*M_PI*x*x*y0)/(x*(1.0 - (2.0*y2)/(x)));
        y[2] = 4*M_PI*x*x*c[0];
        y[3] = (y3)/(x*sqrt(1.0-(2.0*y2)/x));
    }
    //This funciton is not guaranteed to work in all cases. For instance, we have manually 
    //made an exception for x=0, since evaluating at 0 produces infinities and NaNs. 
    //Be sure to declare any exceptions before running, both here and in exceptionHandler(), depending 
    //on the kind of exception desired.  
}


void getInitialCondition (double y[])
{
    //be sure to have these MATCH the equations in diffyQEval
    y[0] = 0.016714611225000002; //Pressure
    y[1] = 0.0; //nu
    y[2] = 0.0; //mass
    y[3] = 0.0; //r-bar
}

void constEval (double y[], double c[])
{
    //Sometimes we want to evaluate constants in the equation that change, but do not have derivative forms.
    //Today, we do that for the total energy density. 
    double c0 = c[0];
    //Make sure to instantiate buffer values so you don't end up changing things when you don't want to!
    //Not necessary here, but left in as a demonstration. 
    c[0] = sqrt(y[0]) + y[0];
    //The total energy density only depends on pressure. 
}


void knownQEval (double x, double y[])
{
    //This function is only used if there are known solutions. 
    //Notably this is not the case for the TOV equations. 
    //If you do put anything here, make SURE it has the same order as the differential equations. 
    //In the case of TOV, that would be Pressure, nu, mass, and r-bar, in that order. 
}


void exceptionHandler (double x, double y[], double c[])
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


int doWeTerminate (double x, double y[], double c[])
{
    //This funciton might be empty. It's only used if the user wants to have a special termination condition.
    //Today we do. We terminate once the pressure hits zero, or goes below it. 
    if (y[0] <= 0.0) {
        return 1;
    } else {
        return 0;
    }
    //return 1 for termination.
}


/*
 * Complicated Example: TOV Solver With Basic Assumptions.
 */
int main() {

double butcher[5][5] = {{0.0,0,0,0,0},{0.5,0.5,0,0,0},{0.5,0.0,0.5,0,0},{1.0,0.0,0.0,1.0,0},{4.0,0.16666666666666666,0.3333333333333333,0.3333333333333333,0.16666666666666666}};
    printf("Beginning ODE Solver \"Odie\" V7...\n");
    
    //SECTION I: Preliminaries
    //Before the program actually starts, variables need to be created
    //and set, as well as the functions chosen. 
    //The system of differential equations can be found declared in diffyQEval().

    double step = 0.00001; //the "step" value.
    double bound = 0.0; //where the boundary/initial condition is. Same for every equation in the system.
    int numberOfEquations = 4; //How many equations are in our system?
    int numberOfConstants = 1; //How many constants do we wish to separately evaluate and report? 
    //If altering the two "numberOf" ints, be careful it doesn't go over the actual number and cause an overflow 
    //in the functions above main()
    const int SIZE = 100000; //How many steps we are going to take?
    bool validate = false; //Set to true if you wish to run a validation test. Only works if solution is already known.
    //Spits out nonsense if no solution is provided.
    //BE WARNED: setting validate to true makes it print out all error data on a second line, the file will have
    //to be read differently.

        //How to get array size: https://stackoverflow.com/questions/37538/how-do-i-determine-the-size-of-my-array-in-c
    size_t methodSize = sizeof(butcher)/sizeof(butcher[0][0]);
    int dimension = sqrt((int)methodSize);
    //We need to know how big our method is, especially if passed one we've never seen before. 
    if (validate == true) {
        printf("Method Order: %i. \nOrder of Error should be near to or larger than Method Order + 1.\n",(int)butcher[dimension-1][0]);
        printf("If not, try a larger step size, roundoff error may be interfering.\n");
    } else {
        printf("Method Order: %i.\n",(int)butcher[dimension-1][0]);
    }
    //If validation is not needed, we don't care about the Order of the Error. 

    double y[numberOfEquations];
    double c[numberOfConstants];
    //These variables temporarily store the values calculated before they are 
    //printed to the output file and forgotten.
    //y is the values of the actual equations. 
    //c is just used to hold any constants we wish to report. 
    //Each array only holds values at one evaluation point, but one for each Equation.

    //This here sets the initial conditions as declared in getInitialCondition()
    getInitialCondition(y); 

    //This evaluates any constants that might be needed for evaluating the actual differnetial equations. 
    constEval(y,c);

    //SECTION II: The Loop

    //prior to beginning the loop, start the timer. 
    double startN = time(NULL);
    //printf("Time: %.90f seconds\n", startN);
    //set start time to current time.  Uncomment to print.

    //also open the file we'll be writing data to. 
    FILE *fp;
    fp = fopen("oCData.txt","w");

    //First though, let's print out our initial data. The print function needs to be adaptable to any size of data. 
    //We can do this with multiple print functions and just not adding the newline character until we're done.
    //We print both to console and to the file for the initial conditions, but later only print to file.
    //First, print the location we are at. 
    printf("INITIAL: Position:,\t%f,\t",bound);
    fprintf(fp, "Position:,\t%f,\t",bound);
    //Second, go through and print the result for every single equation in our system.
    for (int n = 0; n < numberOfEquations; n++) {
        printf("Equation %i:,\t%10.9e,\t",n, y[n]);
        fprintf(fp, "Equation %i:,\t%10.9e,\t",n, y[n]);
    }
    //Third, print out desired constants.     
    for (int n = 0; n < numberOfConstants; n++) {
        printf("Constant %i:,\t%10.9e,\t",n, c[n]);
        fprintf(fp, "Constant %i:,\t%10.9e,\t",n, c[n]);
    }
    //Lastly, the newline character. 
    printf("\n");
    fprintf(fp,"\n");
    //Comma delimiters are printed to the file so it can be converted to .csv with ease. 

    if (validate == true) {
        //In order to keep things neat and regular in the file, print a first line of errors. 
        //Even though by necessity all of them must be zero. 
        fprintf(fp, "Errors:,\t");
            for (int n = 0; n < numberOfEquations; n++) {
                fprintf(fp, "Equation %i:,\t0.0,\t",n);
                fprintf(fp, "Truth:,\t%10.9e,\t",y[n]);
            }
            for (int n = 0; n < numberOfConstants; n++) {
                fprintf(fp, "Constant %i:,\t0.0,\t",n);
                fprintf(fp, "Truth:,\t%10.9e,\t",c[n]);
            }   
            //printf("\n");
            fprintf(fp,"\n");
    }
    
    //This loop fills out all the data.
    //It takes a provided butcher table and executes the method stored within. Any table should work. 

    for (int i = 0; i < SIZE; i++){ 
        //i represents how many steps have been taken. 0 is the initial condition, that is, the variable `bound`. 
        double K[dimension][numberOfEquations];
        //These are the K-values that are required to evaluate RK-like methods. 
        //They will be determined based on the provided butcher table.
        //This is a 2D matrix since each diffyQ has its own set of K-values. 

        //Since we'll be calling K while it's empty, even though there should be no errors due
        //to the way it's set up, let's go ahead and fill it with zeroes.
        for (int j = 0; j<dimension; j++) {
            for (int n = 0; n<numberOfEquations; n++) {
                K[j][n]=0.0;
            }
        } 

        double yInsert[numberOfEquations];
        //We also need an array for the inserted y-values for each equation. 
        //Most applications actually have the different yInsert values be independent, so 
        //if we knew the form of the equation we could simplify the code.
        //However, we need to make sure to always fill everything in case we have a system
        //of the form y'=f(u,y) u'=g(u,y)

        double cInsert[numberOfConstants];
        //Create an array to hold the constants we want.
        //Surprisingly, does not throw an error when the number is zero. Neat.

        for (int j = 1; j < dimension; j++) {
            //Due to the way the Butcher Table is formatted, start our index at 1 and stop at the end. 
            double xInsert = bound+i*step + butcher[j-1][0]*step;
            //x does not change much for different tables, just adjust the "step correction" term.
            //x is the same for every equation too.

            for (int n = 0; n < numberOfEquations; n++) {
                yInsert[n] = y[n];
            } 
            //This is set since y is our actual value, but we will be adjusting yInsert a lot to find the next y. 

            for (int n = 1; n < dimension; n++) {
                //Once again, start at index of 1 rather than 0.
                for (int q = 0; q < numberOfEquations; q++) {
                    yInsert[q] = yInsert[q] + butcher[j-1][n]*K[n][q];
                }
                //Each individual yInsert portion is dependent on one of the K values.
                //K values are initially set to zero even though technically whenever 
                //we would use an undeclared K-value the butcher table would have zero.
                //You know, just in case something goes wrong. 
            }
            
            //Check for any limitations on our results. 
            exceptionHandler(xInsert,yInsert,cInsert);

            //Evaluate the constants. 
            constEval(yInsert,cInsert);

            //Now we actually evaluate the differential equations. 
            diffyQEval(xInsert, yInsert, cInsert);
            //yInsert comes out as evaluated derivatives, and is now in a form we can use. 

            for (int n = 0; n < numberOfEquations; n++) {
                K[j][n] = step*yInsert[n];
                //Fill in the K-values. 
            } 
        }

        //Now that we have all the K-values set, we need to find the actual result in one final loop.
        for (int n = 0; n< numberOfEquations; n++) {
            K[0][n] = y[n]; //The 0th spot in the K-values is reserved for holding the 
            //final value while it's being calculated. 
            for (int j = 1; j < dimension; j++) {
                K[0][n] = K[0][n] + butcher[dimension-1][j]*K[j][n]; 
                //This is where the actual approximation is finally performed. 
            }
            if (validate == true) {
                yInsert[n] = y[n];
                //Before we change our initial ys, we save their data for the validation check now that
                //yInsert is not being used. 
            }
            y[n] = K[0][n]; //Set y to the new estimated value. 
        }

        //After each step is calculated, print results. 
        //However, prior to printing we need to run our exception and constant evaluators one more time. 
        exceptionHandler(bound+i*step,y,c);
        constEval(y,c);
        //Since we've usually been running them on yInsert, the actual y and c values have not generally seen 
        //the restrictions applied here. 

        //Uncomment for live updates.
        //printf(fp, "Position:,\t%f,\t",bound+(i+1)*step);
        fprintf(fp, "Position:,\t%f,\t",bound+(i+1)*step);
        for (int n = 0; n < numberOfEquations; n++) {
            //printf("Equation %i:,\t%10.9e,\t",n, y1[n]);
            fprintf(fp, "Equation %i:,\t%10.9e,\t",n, y[n]);
        }
        for (int n = 0; n < numberOfConstants; n++) {
            //printf("Constant %i:,\t%10.9e,\t",n, c[n]);
            fprintf(fp, "Constant %i:,\t%10.9e,\t",n, c[n]);
        }   
        //printf("\n");
        fprintf(fp,"\n");
                
        //validation: grab the first error for every equation, estimate its order, and continually report errors
        if (validate==true) {
            //We should only be here if we have a truth to compare against. 
            double saveErr1[numberOfEquations];
            double yTruth[numberOfEquations];
            double constErr1[numberOfConstants];
            double cTruth[numberOfConstants];
            //we need errors for both the differential equations and the constants. 
            //We only need one error array for reporting errors, but we will need a second for order checking. 
            exceptionHandler(bound+step*(i+1),yTruth,cTruth);
            knownQEval(bound+step*(i+1),yTruth);
            constEval(yTruth,cTruth);
            //Evaluate what yTruth is at our current point.
            //remember, i+1, not i, since we've already updated the previous point. 

            for (int n = 0; n < numberOfEquations; n++) {
                saveErr1[n] = (yTruth[n] - y[n]);
                //This calculates the errors we have right now. 
            } 
            for (int n = 0; n < numberOfConstants; n++) {
                constErr1[n] = (cTruth[n] - c[n]);
                //This calculates the constant errors we have right now. 
            } 

            if(i == 0.0) { 
                //The following is an algorithm for determining the rate of error 
                //convergence. A bit rudimentary, could be condensed, but is also only
                //called once so not relaly a concern and it is easier to read this way. 
                //Note that this only reports the estimated order for the differential equations, not constants. 
                double saveErr2[numberOfEquations];
                //Need to store a second error. 
                double step2 = step*0.5;
                //It is easier to just use another variable than multiply the step by 0.5 every time. 
                //Below not really commented since it is just a copy of what's above with minor tweaks.
                    for (int j = 0; j<dimension; j++) {
                        for (int n = 0; n<numberOfEquations; n++) {
                            K[j][n]=0.0;
                        }
                    }

                    double yInsertBuffer[numberOfEquations];

                    for (int j = 1; j < dimension; j++) {
                        double xInsert = bound+i*step2 + butcher[j-1][0]*step2;
                        
                        for (int n = 0; n < numberOfEquations; n++) {
                            yInsertBuffer[n] = yInsert[n];
                        } 
                        for (int n = 1; n < dimension; n++) {
                            for (int q = 0; q < numberOfEquations; q++) {
                                yInsertBuffer[q] = yInsertBuffer[q] + butcher[j-1][n]*K[n][q];
                            }
                        }

                        exceptionHandler(xInsert,yInsertBuffer,cInsert);

                        constEval(yInsertBuffer,cInsert);

                        diffyQEval(xInsert, yInsertBuffer, cInsert);
                        for (int n = 0; n < numberOfEquations; n++) {
                            K[j][n] = step2*yInsertBuffer[n];
                        }

                    }
                    for (int n = 0; n< numberOfEquations; n++) {
                        K[0][n] = yInsert[n];
                        for (int j = 1; j < dimension; j++) {
                            K[0][n] = K[0][n] + butcher[dimension-1][j]*K[j][n];
                        }
                    }
                //Now that we've performed the approximation's first step at half the size, we can estimate the order.
                //Create an array to hold the true values. 
                double truthValidate[numberOfEquations];
                //Fill it with the true values. 
                knownQEval(bound+step2,truthValidate);
                //Then print out the estimated order, one individually for each equation.
                for (int n = 0; n < numberOfEquations; n++) {
                    saveErr2[n] = (truthValidate[n] - K[0][n]);
                    printf("Order of Error: %i\t%f\n",n, log2(saveErr1[n]/saveErr2[n]));
                }
                //Note: this will not produce an integer, but with proper data it will be close to an integer
                //and the validation would be performed by rounding. 
                //Results can be larger than the order of the method+1, but should not be much smaller
                //However one can also get errors if the results are too exact, roundoff error can ruin the calcluation. 
                //Using larger step sizes usually removes that. 
            }
        
            //When validating, we add an extra line to the output file for errors.
            //prints the errors alongside the truth values. 
            fprintf(fp, "Errors:,\t");
            for (int n = 0; n < numberOfEquations; n++) {
                fprintf(fp, "Equation %i:,\t%10.9e,\t",n, saveErr1[n]);
                fprintf(fp, "Truth:,\t%10.9e,\t",yTruth[n]);
            }
            for (int n = 0; n < numberOfConstants; n++) {
                fprintf(fp, "Constant %i:,\t%10.9e,\t",n, constErr1[n]);
                fprintf(fp, "Truth:,\t%10.9e,\t",cTruth[n]);
            }   
            //printf("\n");
            fprintf(fp,"\n");

            //Note that error printed is not relative error, but simple difference error. 
        }

        //And the very last thing we do in the loop is ask if we terminate it. 
        if (doWeTerminate(bound+i*step, y, c) == 1) {
            i = SIZE;
        }
    }

    //SECTION III: Analysis
    //Minor post-processing goes here. 
    //Anything advanced will need to be done in the jupyter notebook that encloses this. 

    // basic reference: https://www.tutorialspoint.com/cprogramming/c_file_io.htm
    // used to be a file converter here, now there isn't, we just close the file. 
    fclose(fp);

    //TIMER
    double endN = time(NULL);
    //loop is complete, how long did it take?
    printf("Time Elapsed: %f seconds\n", endN-startN);
    //Only calculates to the closest second, for some reason.

    printf("ODE Solver \"Odie\" V7 Shutting Down...\n");
    return 0;

// - GM, master of dogs.
    }
