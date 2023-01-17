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
//This fifth version seeks to allow the program to solve a system of ODEs, starting with a very simple one.
//User functionality will (eventually) involve input of own equations, but we don't have that here yet, and it may be a while.

//Heavily influenced by Numerical Mathematics and Computing 6E by Cheney and Kincaid.

//Outside the program, we substantiate the differential equation itself.

double diffyQEval (double x, double y[], int i)
{
    switch(i){
        case 0: {
            return -((1.0+y[0])*( (2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0] ))/(x*2.0*(1.0 - (2.0*y[2])/(x))); 
            break;
        }
        case 1: {
            return ((2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0])/(x*(1.0 - (2.0*y[2])/(x)));
            break;
        }
        case 2: {
            return 4*M_PI*x*x;
            break;
        }
        case 3: {
            return (y[3])/(x*sqrt(1.0-(2.0*y[2])/x));
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

double getInitialCondition (int n)
{
    //be sure to have these MATCH the equations in diffyQEval
    switch(n){
        case 0: {
            return 0.016714611225000002; 
            break;
        }
        case 1: {
            return 0.0;
            break;
        }
        case 2: {
            return 6.115617413201306e-16;
            break;
        }
        case 3: {
            return 6.998108874980957e-06;
            break;
        }
        // case...
        //add more cases for higher order systems. 
        default: {
            //Any number could be a result, theoretically, so there isn't an easy error value to put here.
        }
    }
}

//Remember when adjusting these to adjust the boundary value bValue in main() as well. 

int main()
{
    printf("Beginning ODE Solver \"Odie\" V4...\n");
    
    //SECTION I: Preliminaries
    //Before the program actually starts, variables need to be created
    //and set, as well as the function itself chosen. 
    //The diffyQ itself can be found declared in diffyQEval().

    //Butcher Table: for now we define our method table here. 
    //When run through the notebook this section is absent as it fills it itself. 
    //Uncomment the method you wish to use. 
    //double butcher[4][4] = {{0.0,0.0,0.0,0.0},{1.0,1.0,0.0,0.0},{0.5,0.25,0.25,0.0},{3,1.0/6.0,1.0/6.0,2.0/3.0}};
    //This is the SSPRK3 method, chosen since it has a simple array but with less zeroes than other options. 
    //Note that the "3" in the last row is the order, that slot is always empty on a butcher table. 
    //If you wish for a different method, comment out the one above and initate one of the below:

    //double butcher[2][2] = {{0.0,0.0},{1,1.0}};
    //This is Euler's Method, good for test cases since it's easy to look at. 

    //double butcher[3][3] = {{0.0,0.0,0.0},{1.0,1.0,0.0},{2,0.5,0.5}};
    //RK2

    //double butcher[5][5] = {{0.0,0.0,0.0,0.0,0.0},{0.5,0.5,0.0,0.0,0.0},{0.5,0.0,0.5,0.0,0.0},{1.0,0.0,0.0,1.0,0.0},{4,1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0}};
    //RK4

    double butcher[7][7] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0},
    {0.2,0.2,0.0,0.0,0.0,0.0,0.0},
    {0.3,3.0/40.0,9.0/40.0,0.0,0.0,0.0,0.0},
    {0.6,0.3,-9.0/10.0,1.2,0.0,0.0,0.0},
    {1.0,-11.0/54.0,2.5,-70.0/27.0,35.0/27.0,0.0,0.0},
    {7.0/8.0,1631.0/55296.0,175.0/512.0,575.0/13824.0,44275.0/110592.0,253.0/4096.0,0.0},
    {5,37.0/378.0,0.0,250.0/621.0,125.0/594.0,0.0,512.0/1771.0}};
    //RK5 (Cash-Karp version)

    //How to get array size: https://stackoverflow.com/questions/37538/how-do-i-determine-the-size-of-my-array-in-c
    size_t methodSize = sizeof(butcher)/sizeof(butcher[0][0]);
    int dimension = sqrt((int)methodSize);
    printf("Method Order: %i. \nOrder of Error should be near Method Order + 1.\n",(int)butcher[dimension-1][0]);
    printf("If not, try a larger step size, roundoff error may be interfering.\n");

    double step = 0.001; //the "step" value.
    double bound = 0.00000001; //where the boundary/initial condition is. Same for every equation in the system.
    int numberOfEquations = 4; //How many equations are in our system?
    const int SIZE = 2000; //How many steps we are going to take?
    bool validate = false; //set to true if you wish to run a validation test.
    //Attempts to find the order of the method used. 

    double bValue[numberOfEquations]; 

    //This here sets the initial conditions as declared in getInitialCondition()
    for (int n = 0; n < numberOfEquations; n++) {
        bValue[n] = getInitialCondition(n);
    }


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

    //First though, let's print out our data. The print function needs to be adaptable to any size of data. 
    //We can do this with multiple print functions and just not adding the newline character until we're done.
    //Commented out version is for if we want to compare with a known function.
    //printf("INITIAL: Position:\t%f\tTruth:\t%10.9e\tCalculated:\t%10.9e\tError:\t%10.9e\t\n",bound,yTruth1, y1[0], yError);  
    printf("INITIAL: Position:,\t%f,\t",bound);
    fprintf(fp, "Position:,\t%f,\t",bound);
    for (int n = 0; n < numberOfEquations; n++) {
        printf("Equation %i:,\t%10.9e,\t",n, y1[n]);
        fprintf(fp, "Equation %i:,\t%10.9e,\t",n, y1[n]);
    }
    printf("\n");
    fprintf(fp,"\n");

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
        //uncomment printf sections if you want live updates. 
        //printf("Position:,\t%f,\t",bound);
        fprintf(fp, "Position:,\t%f,\t",bound+(i+1)*step);
        for (int n = 0; n < numberOfEquations; n++) {
            //printf("Equation %i:,\t%10.9e,\t",n, y1[n]);
            fprintf(fp, "Equation %i:,\t%10.9e,\t",n, y1[n]);
        }
        //printf("\n");
        fprintf(fp,"\n");
                
        //validation: grab the first nonzero error, calculate its order.
        //Currently broken. 
        if(validate==true && i == 0.0) { 
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
}

// - GM, master of dogs.