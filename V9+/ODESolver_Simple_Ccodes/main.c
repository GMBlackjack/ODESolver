#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdbool.h"
#include "time.h"

void diffyQEval (double x, double y[], double c[])
{
    //Efficient(ish) Assignment Method

        double y0 = y[0];
        double y1 = y[1];
        y[0] = y1;
        y[1] = y0 + x;

    //This is the differential equation system itself. 
    //We have a very simple y'' = y+x situation here, split up into
    // y[0]' = y[1]
    // y[1]' = y[0]+x
    //Naturally other equaitons can be put in, but be sure to change the numberOfEquations value!
    //Note: not guaranteed to work for functions that are not well-behaved. 
}


void getInitialCondition (double y[])
{
    //be sure to have these MATCH the equations in diffyQEval
    y[0] = 2.0;
    y[1] = -1.0;
}


void constEval (double y[], double c[])
{
    //none. 
}


void knownQEval (double x, double y[])
{
    y[0] = exp(x) + exp(-x) - x;
    y[1] = exp(x) - exp(-x) - 1;
    //This function is only used if there are known solutions. 
    //Do note that this would change with different boundary conditions. 
}



void exceptionHandler (double x, double y[], double c[])
{
    //This funciton might be empty. It's only used if the user wants to hard code some limitations 
    //On some varaibles.
    //Good for avoding divide by zero errors, or going negative in a square root. 
}

int doWeTerminate (double x, double y[], double c[])
{
    //This funciton might be empty. It's only used if the user wants to have a special termination condition.
    return 0;
    //return 1 for termination.
}

/*
 * Simple Example: u''=u+x Solver
 */
int main() {

double butcher[15][14] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{1.0/18.0,1.0/18.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{1.0/12.0,1.0/48.0,1.0/16.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{1.0/8.0,1.0/32.0,0.0,3.0/32.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{5.0/16.0,5.0/16.0,0.0,-75.0/64.0,75.0/64.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{3.0/8.0,3.0/80.0,0.0,0.0,3.0/16.0,3.0/20.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{59.0/400.0,29443841.0/614563906.0,0.0,0.0,77736538.0/692538347.0,-28693883.0/1125000000.0,23124283.0/1800000000.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{93.0/200.0,16016141.0/946692911.0,0.0,0.0,61564180.0/158732637.0,22789713.0/633445777.0,545815736.0/2771057229.0,-180193667.0/1043307555.0,0.0,0.0,0.0,0.0,0.0,0.0},{5490023248.0/9719169821.0,39632708.0/573591083.0,0.0,0.0,-433636366.0/683701615.0,-421739975.0/2616292301.0,100302831.0/723423059.0,790204164.0/839813087.0,800635310.0/3783071287.0,0.0,0.0,0.0,0.0,0.0},{13.0/20.0,246121993.0/1340847787.0,0.0,0.0,-37695042795.0/15268766246.0,-309121744.0/1061227803.0,-12992083.0/490766935.0,6005943493.0/2108947869.0,393006217.0/1396673457.0,123872331.0/1001029789.0,0.0,0.0,0.0,0.0},{1201146811.0/1299019798.0,-1028468189.0/846180014.0,0.0,0.0,8478235783.0/508512852.0,1311729495.0/1432422823.0,-10304129995.0/1701304382.0,-48777925059.0/3047939560.0,15336726248.0/1032824649.0,-45442868181.0/3398467696.0,3065993473.0/597172653.0,0.0,0.0,0.0},{1.0,185892177.0/718116043.0,0.0,0.0,-3185094517.0/667107341.0,-477755414.0/1098053517.0,-703635378.0/230739211.0,5731566787.0/1027545527.0,5232866602.0/850066563.0,-4093664535.0/808688257.0,3962137247.0/1805957418.0,65686358.0/487910083.0,0.0,0.0},{1.0,403863854.0/491063109.0,0.0,0.0,-5068492393.0/434740067.0,-411421997.0/543043805.0,652783627.0/914296604.0,11173962825.0/925320556.0,-13158990841.0/6184727034.0,3936647629.0/1978049680.0,-160528059.0/685178525.0,248638103.0/1413531060.0,0.0,0.0},{8.0,14005451.0/335480064.0,0.0,0.0,0.0,0.0,-59238493.0/1068277825.0,181606767.0/758867731.0,561292985.0/797845732.0,-1041891430.0/1371343529.0,760417239.0/1151165299.0,118820643.0/751138087.0,-528747749.0/2220607170.0,1.0/4.0},{8.0,13451932.0/455176623.0,0.0,0.0,0.0,0.0,-808719846.0/976000145.0,1757004468.0/5645159321.0,656045339.0/265891186.0,-3867574721.0/1518517206.0,465885868.0/322736535.0,53011238.0/667516719.0,2.0/45.0,0.0}};
    printf("Beginning ODE Solver \"Odie\" V7...\n");
    
    //SECTION I: Preliminaries
    //Before the program actually starts, variables need to be created
    //and set, as well as the functions chosen. 
    //The system of differential equations can be found declared in diffyQEval().

    double step = 0.05; //the "step" value.
    double bound = 0.0; //where the boundary/initial condition is. Same for every equation in the system.
    int numberOfEquations = 2; //How many equations are in our system?
    int numberOfConstants = 0; //How many constants do we wish to separately evaluate and report? 
    //If altering the two "numberOf" ints, be careful it doesn't go over the actual number and cause an overflow 
    //in the functions above main()
    const int SIZE = 20; //How many steps we are going to take?
    bool validate = true; //Set to true if you wish to run a validation test. Only works if solution is already known.
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
    fp = fopen("oSData.txt","w");

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
