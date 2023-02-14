//#include "nrpy_odiegm.h"
#include "butcher.c"
//TODO:
//GSL Functionality
//Jupyter notebook Adjusents. 

//Note: math.h requries the "-lm" arg be added at the END of tasks.json's arguments.
//https://askubuntu.com/questions/332884/how-to-compile-a-c-program-that-uses-math-h

//ODE Solver "Odie"
//By G. M. Steward
//The main goal of this project is to solve Ordinary Differential Equation Systems
//in complete generality. However, we will take a while to get there.
//This ninth version seeks to make this code functional as a drop-in replacement for GSL's solver. 

//Heavily influenced by Numerical Mathematics and Computing 6E by Cheney and Kincaid.
//and GSL's ODE Solver, especially the method for adaptive time step and high-level funcitonality. 

//Outside the program, we substantiate the differential equation itself as well as its various limitations. 

//We have to define this outside so we can pass it to various functions.
struct constantParameters { 
    int dimension; //number that says how many we have. 
    double rho;
    //add more as necessary.  
};

void exceptionHandler (double x, double y[], struct constantParameters *params)
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

int doWeTerminate (double x, double y[], struct constantParameters *params)
{
    //This funciton might be empty. It's only used if the user wants to have a special termination condition.
    //Today we do. We terminate once the pressure hits zero, or goes below it. 
    if (y[0] < 1e-16) {
        return 1;
    } else {
        return 0;
    }
    //return 1 for termination.
}

void constEval (double x, const double y[], struct constantParameters *params)
{
    //Sometimes we want to evaluate constants in the equation that change, but do not have derivative forms.
    //Today, we do that for the total energy density. 
    params->rho = sqrt(y[0]) + y[0];
    //The total energy density only depends on pressure. 
}

int diffyQEval (double x, const double y[], double dydx[], void *params)
{
    //GSL-adapted evaluation function. 
    //It is possible to do this with one array, but GSL expects two. 

    //dereference the struct
    double rho = (*(struct constantParameters*)params).rho;
    //WHY oh WHY GSL do you demand we use a VOID POINTER to the struct...
    //https://stackoverflow.com/questions/51052314/access-variables-in-struct-from-void-pointer

    //This if statement is an example of a special condition, in this case at x=0 we have a divide by zero problem. 
    //In this case we manually know what the derivatives should be.
    //Alternatively, we could define piecewise equations this way. 
    if(x == 0) {
        dydx[0] = 0; 
        dydx[1] = 0;
        dydx[2] = 0;
        dydx[3] = 1;
    }
    else {
        dydx[0] = -((rho+y[0])*( (2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0] ))/(x*2.0*(1.0 - (2.0*y[2])/(x)));
        dydx[1] =  ((2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0])/(x*(1.0 - (2.0*y[2])/(x)));
        dydx[2] = 4*M_PI*x*x*rho;
        dydx[3] = (y[3])/(x*sqrt(1.0-(2.0*y[2])/x));
    }
    //This funciton is not guaranteed to work in all cases. For instance, we have manually 
    //made an exception for xdouble butcher[3][3] = {{0.0,0,0},{1.0,1.0,0},{2.0,0.5,0.5}};=0, since evaluating at 0 produces infinities and NaNs. 
    //Be sure to declare any exceptions before running, both here and in exceptionHandler(), depending 
    //on the kind of exception desired.  

    return 1;
}

//This is the function to evaluate the known solution. Must be set manually.
void knownQEval (double x, double y[])
{
    //y[0] = ...
    //y[1] = ...
    //This function is only used if there are known solutions. 
    //Notably this is not the case for the TOV equations. 
    //If you do put anything here, make SURE it has the same order as the differential equations. 
    //In the case of TOV, that would be Pressure, nu, mass, and r-bar, in that order. 
}

void getInitialCondition (double y[])
{
    //be sure to have these MATCH the equations in diffyQEval
    y[0] = 0.016714611225000002; //Pressure, can be calcualated from central baryon density. 
    y[1] = 0.0; //nu
    y[2] = 0.0; //mass
    y[3] = 0.0; //r-bar
}

void assignConstants (double c[], struct constantParameters *params)
{
    //this really should be handled by the computer automatically, we just don't know how.
    c[0] = params->rho;
    //add more as required. 
}

//Remember when adjusting these to adjust the necessary parameters in main() as well:
//step, bound, numberOfEquations, numberOfConstants, SIZE, and validate. 

int main()
{
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
    const int SIZE = 100000; //How many steps we are going to take? This is the default termination condition. 
    bool validate = false; //Set to true if you wish to run a validation test. Only works if solution is already known.
    //Spits out nonsense if no solution is provided.
    //BE WARNED: setting validate to true makes it print out all error data on a second line, the file will have
    //to be read differently. 
    bool reportErrorEstimates = false;
    //enable to have every other line in "oData.txt" report the estimated erorr values.
    bool reportErrorActual = false;
    //enable to have the actual errors reported. Produces junk data if there is no declared known function. 

    //ERROR PARAMETERS: Use these to set limits on the erorr. 
    double absoluteErrorLimit = 1e-10; //how big do we let the absolute error be?
    double relativeErrorLimit = 1e-10; //how big do we let the relative error be?
    double ayErrorScaler = 1.0; //For giving extra weight to the functions themselves.
    double adyErrorScaler = 1.0;  //For giving extra weight to the derivatives. 
    //Note: the function for determining error creates a combined tester, it doesn't do either/or absolute or relative.
    //So its' best to either cap absolute or relative error, not both. 
    //Constants should be set to 1 unless you know what you're dong, which I don't.

    //adaptive timestep constraints
    //We allow for some truly adaptable control to how the adaptive timestep works. All the paramaters can be adjusted. 
    double scaleFactor = 0.9; //a general scaling factor that keeps the step size from "freaking out." 0.9 by defaut. 
    double errorSafety = 4.0/15.0; //A mysterious factor GSL uses to try to estimate the error. 
    // 4.0/15.0 by default. Not sure where it comes from. 
    double errorUpperTolerance = 1.1; //if the error ratio is larger than this, lower the step size. 1.1 (10% over) by default.
    double errorLowerTolerance = 0.5; //if the error ratio is lower than this, raise the step size. 0.5 (50% under) by default. 
    double maxStepAdjustment = 5.0; //the largest single adjustment that can be made to the step size. 5 default. 
    double minStepAdjustment = 0.2; //the minimum single adjustment that can be made to the step size. 1/5 default. 
    double absoluteMaxStep = 0.1; //An absolute cap on the step size, 0.1 by default. 
    double absoluteMinStep = 1e-10; //An absolute floor on the step size, 1e-10 by default. 

    //double butcher[5][5] = {{0.0,0,0,0,0},{0.5,0.5,0,0,0},{0.5,0.0,0.5,0,0},{1.0,0.0,0.0,1.0,0},{4.0,0.16666666666666666,0.3333333333333333,0.3333333333333333,0.16666666666666666}};

    //We need to define a struct that can hold all possible constants. 
    struct constantParameters cp; 

    cp.dimension = numberOfConstants;
    //we'll set the actual density later. 

    nrpy_odiegm_system system = {diffyQEval,NULL,numberOfEquations,&cp};
    //This is the system of equations we solve.
    //Second null should point to constants. 

    //Now we set up the method. This is definitely a very roundabout way to do it,
    //And for now all the pointer nonsense is immediately undone, but here it is. 
    //All the butcher tables themselves are defined in butcher.c. 
    //We just need to create an object that gets them.
    const nrpy_odiegm_step_type * stepType;
    stepType = nrpy_odiegm_step_AB;
    //Here is where the method is actually set, by specific name since that's what GSL does. 

    int methodType = 1;
    if (stepType->rows == stepType->columns) {
        methodType = 0; //aka, normal RK-type method. 
    } //no need for an else, we set it to 1 earlier to represent Adaptive methods. 
    //This integer is actually very useful as it can help us keep track of index changes between
    //the two types of methods! 

    //Technically this sets the type. What we do now is basically undo that pointer nonsense. 
    int rows = stepType->rows;
    int columns = stepType->columns;
    //Since we know the dimension, we can now fill an actual butcher array. 
    if (rows == 1) {
        //We're using the Adams method, special table construction required. 
        methodType = 2;
        rows = 5;
        columns = 5;
        //in the future this number can be specified, but not today. 
    }
    double butcher[rows][columns];
    int counter = 0;
    if (methodType == 2) {
        //Adams methods have a very different implementation. However it is much simpler once the table is constructed!
        //For now, let's just explicitly grab the 5th order method. 
        //Ironically the best way to do that is manually declare an array and then feed it back into butcher. 
        double adamsArray[5][5] = {{1901.0/720.0,-2774.0/720.0,2616.0/720.0,-1274.0/720.0,251.0/720.0},{0,55.0/24.0,-59.0/24.0,37.0/24.0,-9.0/24.0},{0,0,23.0/12.0,-16.0/12.0,5.0/12.0},{0,0,0,3.0/2.0,-1.0/2.0},{5,0,0,0,1.0}};
        for (int n=0; n < rows; n++) {
            for (int m = 0; m < columns; m++) {
                butcher[n][m] = adamsArray[n][m];
                printf("%10.9e  ", butcher[n][m]);
            }
            printf("\n");
        }
    } else {
        for (int i=0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                butcher[i][j] = *((double *)(*stepType).butcher+counter);
                counter++;
            }
        }
    }

    //How to get array size no longer needed.

    if (validate == true) {
        printf("Method Order: %i. \nOrder of Error should be near to or larger than Method Order + 1.\n",(int)butcher[rows-1][0]);
        printf("If not, try a larger step size, roundoff error may be interfering.\n");
    } else {
        printf("Method Order: %i.\n",(int)butcher[rows-1][0]);
    }
    //If validation is not needed, we don't care about the Order of the Error. 

    double y[numberOfEquations];
    //These variables temporarily store the values calculated before they are 
    //printed to the output file and forgotten.
    //y is the values of the actual equations. 
    //c is just used to hold any constants we wish to report. 
    //Each array only holds values at one evaluation point, but one for each Equation.

    double c[numberOfConstants];
    //You'd think that, since we have the constants in a struct, we can avoid declaring this.
    //No. Not as far as we can tell, anyway. Structs are a pain to iterate through,
    //and we can't know what form the user is going to hand us the struct in. 

    double currentPosition = bound;
    //since we adjust the step size it is not possible to algorithmically determine our position. 
    //Thus this variable is required.

    //This here sets the initial conditions as declared in getInitialCondition()
    getInitialCondition(y); 

    //This evaluates any constants that might be needed for evaluating the actual differnetial equations. 
    constEval(currentPosition, y, &cp);

    //SECTION II: The Loop

    //Open the file we'll be writing data to. 
    FILE *fp;
    fp = fopen("oData.txt","w");

    //First though, let's print out our initial data. The print function needs to be adaptable to any size of data. 
    //We can do this with multiple print functions and just not adding the newline character until we're done.
    //We print both to console and to the file for the initial conditions, but later only print to file.
    //First, print the location we are at. 
    printf("INITIAL: Position:,\t%f,\t",bound);
    fprintf(fp, "Position:,\t%f,\t",bound);
    //Second, go through and print the result for every single equation in our system.
    for (int n = 0; n < numberOfEquations; n++) {
        printf("Equation %i:,\t%15.14e,\t",n, y[n]);
        fprintf(fp, "Equation %i:,\t%15.14e,\t",n, y[n]);
    }
    //Third, print out desired constants.
    assignConstants(c,&cp);    
    for (int n = 0; n < numberOfConstants; n++) {
        printf("Constant %i:,\t%15.14e,\t",n, c[n]);
        fprintf(fp, "Constant %i:,\t%15.14e,\t",n, c[n]);
    }
    //Lastly, the newline character. 
    printf("\n");
    fprintf(fp,"\n");
    //Comma delimiters are printed to the file so it can be converted to .csv with ease. 

    if (reportErrorEstimates == true) {
        //In order to keep things neat and regular in the file, print a first line of errors. 
        //Even though by necessity all of them must be zero. 
        fprintf(fp, "Errors Estimates:,\t");
        for (int n = 0; n < numberOfEquations; n++) {
            fprintf(fp, "Equation %i:,\t0.0,\t",n);
        }
        for (int n = 0; n < numberOfConstants; n++) {
            fprintf(fp, "Constant %i:,\t0.0,\t",n);
        }   
        fprintf(fp,"\n");
    }
    
    if (reportErrorActual == true) {
        //In order to keep things neat and regular in the file, print a first line of errors. 
        //Even though by necessity all of them must be zero. 
        fprintf(fp, "Errors:,\t");
        for (int n = 0; n < numberOfEquations; n++) {
            fprintf(fp, "Equation %i:,\t0.0,\t",n);
            fprintf(fp, "Truth:,\t%15.14e,\t",y[n]);
        }
        for (int n = 0; n < numberOfConstants; n++) {
            fprintf(fp, "Constant %i:,\t0.0,\t",n);
            fprintf(fp, "Truth:,\t%15.14e,\t",c[n]);
        }   
        //printf("\n");
        fprintf(fp,"\n");
    }

    //This loop fills out all the data.
    //It takes a provided butcher table and executes the method stored within. Any table should work. 
    if (methodType != 2) {
        //If we're not doing Adams-Bashforth, do the "normal" loop for RK-like methods. 
        for (int i = 0; i < SIZE; i++){ 
            //i represents how many steps have been taken. 0 is the initial condition, that is, the variable `bound`.

            //To use adaptive time-step, we need to store data at different step values:
            double yBigStep[numberOfEquations];
            double ySmolSteps[numberOfEquations];
            //Adams method will consider itself ySmolSteps since RK methods assume that is the "answer."

            double yBigStepHalf[numberOfEquations];
            double ySmolStepsHalf[numberOfEquations];
            //These are here for validation only. We would rather not declare them at all, but they have to be declared outside the loop
            //So the information within can be stored through the iterations. 

            struct constantParameters cpBigStep; 

            cpBigStep.dimension = numberOfConstants;
            cpBigStep.rho = cp.rho;

            struct constantParameters cpSmolSteps; 

            cpSmolSteps.dimension = numberOfConstants;
            cpSmolSteps.rho = cp.rho;
            //Find a way to genrealize this assignment to any number of items in the struct. 

            //One could argue that since the small steps will become our result we shouldn't declare it, however we are actually
            //NOT going to assign them to the actual answer y until we compare and run the adaptive
            //time-step algorithm. We might throw out all the data and need to run it again! 
            double errorEstimate[numberOfEquations+numberOfConstants];
            //even if we aren't limiting the constants, we can still report their error. 
            
            double originalStep = step;
            //We need to be able to refer to the original step so we can see if we're adjusting it too much at once. 

            //We rather explicitly do not actually take any steps until we confirm the error is below what we want.
            bool errorSatisfactory = false;
            bool underError = false;
            bool overError = false;
            //It's important to declare these outside the errorSatisfactory loop since to update the stepper we need to know
            //exactly what kind of step change we just did. 
            while (errorSatisfactory == false) {
                
                //All of the bellow values start off thinking they are the values from the previous step or initial conditions. 
                //We must reset them every time we return here.  
                for (int n = 0; n < numberOfEquations; n++) {
                    yBigStep[n] = y[n];
                    ySmolSteps[n] = y[n];
                } 
                for (int iteration = 1; iteration < 4; iteration++) {
                    //So, we want to use Adaptive Timestep methodology. This will involve evaluating each step three times, 
                    //In order to compare the evolution of two different step sizes and get an error estimate. 
                    //Iteration 1 performs a normal step. 
                    //Iteration 2 perofrms a half step.
                    //Iteration 3 performs another half step after the previous one. 
                    //Naturally the half-step results are reported as truth, but we get an error estimate from the difference
                    //between the two values. 

                    //For adaptive methods we only go through iteration 1 and 2

                    if (i == 0 && iteration == 1 && validate == false && methodType == 0) {
                        //don't take unecessary steps, if we are on the first step and have no need for the large step, ignore it.
                        //Since we always want the first step to go through don't bother calculating things we don't need. 
                        iteration = 2;
                        //This doesn't actually apply to adaptive methods since we cheat and do it in one iteration. 
                    }

                    double scale = 1.0;
                    //this is the number we use to scale. It's either 1 or 1/2, depending on what size step we want. 
                    int shift = 0.0;
                    //this is the number we set if we want to shift where we are evaluating from. 
                    if (iteration == 1.0) {
                        //scale remains 1
                        //shift remains 0
                    } else if (iteration == 2.0) {
                        scale = 0.5; //using half-steps.
                        // shfit remains 0
                    } else {
                        scale = 0.5; //using half-steps.
                        shift = 1; 
                    }
                    //Every time it's needed, we multiply the step by the scale. 

                    double K[rows-methodType][numberOfEquations];
                    //These are the K-values that are required to evaluate RK-like methods. 
                    //They will be determined based on the provided butcher table.
                    //This is a 2D matrix since each diffyQ has its own set of K-values. 
                    //Note that we subtract the method type: K doesn't need to be larger if 
                    //we are using an adaptive butcher table. 

                    //Since we'll be calling K while it's empty, even though there should be no errors due
                    //to the way it's set up, let's go ahead and fill it with zeroes.
                    for (int j = 0; j<rows-methodType; j++) {
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

                    double dyOut[numberOfEquations];
                    //GSL demands that we use two separate arrays for y and y', so here's y'. 

                    double cInsert[numberOfConstants];

                    struct constantParameters cpInsert; 

                    cpInsert.dimension = numberOfConstants;
                    cpInsert.rho = cp.rho;
                    //Create an array to hold the constants we want.
                    //Find way to generalize. 

                    for (int j = 1; j < rows-methodType; j++) {
                        //Due to the way the Butcher Table is formatted, start our index at 1 and stop at the end. 
                        double xInsert = currentPosition+shift*step*scale + butcher[j-1][0]*step*scale;
                        //xInsert does not change much for different tables, just adjust the "step correction" term.
                        //xInsert is the same for every equation too.

                        for (int n = 0; n < numberOfEquations; n++) {
                            yInsert[n] = ySmolSteps[n];
                        } 
                        //Note that we are setting our buffer value, yInsert, to ySmolSteps. 
                        //This is because ySmolSteps is y at first, but we will need to evolve it forward
                        //two steps, so on the second small step this will be different. 

                        for (int n = 1; n < columns; n++) {
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
                        exceptionHandler(xInsert,yInsert,&cpInsert);

                        //Evaluate the constants. 
                        constEval(xInsert, yInsert,&cpInsert);

                        //Now we actually evaluate the differential equations.
                        system.function(xInsert, yInsert, dyOut, &cpInsert);
                        //yInsert comes out as evaluated derivatives, and is now in a form we can use. 

                        for (int n = 0; n < numberOfEquations; n++) {
                            K[j][n] = step*scale*dyOut[n];
                            //Fill in the K-values. 
                        } 
                    }

                    //Now that we have all the K-values set, we need to find the actual result in one final loop.
                    for (int n = 0; n< numberOfEquations; n++) {
                        K[0][n] = ySmolSteps[n]; //The 0th spot in the K-values is reserved for holding the 
                        //final value while it's being calculated. 
                        for (int j = 1; j < columns; j++) {
                            K[0][n] = K[0][n] + butcher[rows-1-methodType][j]*K[j][n]; 
                            //This is where the actual approximation is finally performed. 
                        }
                        ySmolSteps[n] = K[0][n]; //Set ySmol to the new estimated value. 
                    }
                    //Note that we specifically set ySmol to the value, not anything else. 
                    //This is because we wish to avoid abusing if statements and would like to do the following only once.
                    exceptionHandler(currentPosition+(1.0+shift)*step*scale,ySmolSteps,&cpSmolSteps);
                    constEval(currentPosition+(1.0+shift)*step*scale,ySmolSteps,&cpSmolSteps); 
                    //Check for exceptions and evaluate constants. 

                    if (iteration == 1 || (i == 0 && validate == true && methodType == 1)) {
                        for (int n = 0; n<numberOfEquations; n++) {
                            yBigStep[n] = ySmolSteps[n];
                            ySmolSteps[n] = y[n];

                            //we still need to reset the value for SmolSteps on the first iteration
                            //no matter the type of method. 
                        }
                    }
                    //This only runs on the first iteration, setting the big step to the right value
                    //and resetting the small steps for when we actually use it. 
                    //This odd structure exists purely for efficiency. 
                    

                    //If we are in an adaptive method situation, use that method and exit the iterations loop.
                    if (methodType == 1) {
                        for (int n = 0; n< numberOfEquations; n++) {
                            K[0][n] = ySmolSteps[n]; //The 0th spot in the K-values is reserved for holding the 
                            //final value while it's being calculated. 
                            for (int j = 1; j < columns; j++) {
                                K[0][n] = K[0][n] + butcher[rows-1][j]*K[j][n]; 
                                //This is where the actual approximation is finally performed. 
                                //This time we use the bottom row, not the second to bottom row (for adaptive methods)
                            }
                            ySmolSteps[n] = K[0][n]; //Set ySmol to the new estimated value. 
                        }
                        if (validate == true && i == 0 && iteration ==1) {
                            //do nothing fancy. 
                        } else {
                            iteration = 4; //break out now, we don't need to go any further. 
                            //We avoid this on the first iteration only if we want to validate. 
                        }
                    }

                    //printf("test %i %10.9e %10.9e\n", i, yBigStep[0], ySmolSteps[0]);
                    if (validate == true && i == 0) {
                        if (methodType == 0 && iteration == 2) {
                            //The normal validation for the normal RK methods. 

                            //Now that we've performed the approximation's first step at half the size, we can estimate the order.

                            //Create 2 arrays to hold the true values. 
                            double truthValidateBig[numberOfEquations];
                            double truthValidateSmol[numberOfEquations];
                            //Fill it with the true values. 
                            knownQEval(bound+step,truthValidateBig);
                            knownQEval(bound+step*0.5,truthValidateSmol);
                            //Then from this calculate the estimated errors.

                            for (int n = 0; n < numberOfEquations; n++) {
                                truthValidateBig[n] = (truthValidateBig[n] - yBigStep[n]);
                                truthValidateSmol[n] = (truthValidateSmol[n] - ySmolSteps[n]);
                                //Now the validation steps contain their own errors, we can compare them.
                                printf("Order of Error: %i\t%f\n",n, log2(truthValidateBig[n]/truthValidateSmol[n]));
                                //print out the estimated error. 
                            }
                            //Note: this will not produce an integer, but with proper data it will be close to an integer
                            //and the validation would be performed by rounding. 
                            //However one can also get errors if the results are too exact, roundoff error can ruin the calcluation. 
                            //Using larger step sizes usually removes that. 
                        } if (methodType == 1) {
                            if (iteration == 1) {
                                //validation has to be run differently for inherently adaptive methods, since we don't calculate half-steps naturally. 
                                for (int n = 0; n < numberOfEquations; n++) {
                                    yBigStepHalf[n] = yBigStep[n];
                                    ySmolStepsHalf[n] = ySmolSteps[n];
                                    //store the previous values for next loop.
                                    ySmolSteps[n] = y[n];
                                    //This also needs to be cleared again, for once we do not want it to carry over. 
                                }
                            } else {
                                //This should only trigger on iteration == 2. 
                                //The normal validation for the inherently adaptive RK methods. 

                                //Now that we've performed the approximation's first step at half the size, we can estimate the order
                                //for both internal calculations

                                //Create 2 arrays to hold the true values. 
                                double truthValidateBig[numberOfEquations];
                                double truthValidateSmol[numberOfEquations];
                                //Fill it with the true values. 
                                knownQEval(bound+step,truthValidateBig);
                                knownQEval(bound+step*0.5,truthValidateSmol);
                                
                                //Then from this calculate the estimated errors.
                                for (int n = 0; n < numberOfEquations; n++) {
                                    truthValidateBig[n] = (truthValidateBig[n] - yBigStepHalf[n]);
                                    truthValidateSmol[n] = (truthValidateSmol[n] - yBigStep[n]);
                                    //Now the validation steps contain their own errors, we can compare them.
                                    printf("Order of ErrorA: %i\t%f\n",n, log2(truthValidateBig[n]/truthValidateSmol[n]));
                                    //print out the estimated error. 
                                }
                                knownQEval(bound+step,truthValidateBig);
                                knownQEval(bound+step*0.5,truthValidateSmol);

                                for (int n = 0; n < numberOfEquations; n++) {
                                    truthValidateBig[n] = (truthValidateBig[n] - ySmolStepsHalf[n]);
                                    truthValidateSmol[n] = (truthValidateSmol[n] - ySmolSteps[n]);
                                    //Now the validation steps contain their own errors, we can compare them.
                                    printf("Order of ErrorB: %i\t%f\n",n, log2(sqrt(truthValidateBig[n]*truthValidateBig[n])/sqrt(truthValidateSmol[n]*truthValidateSmol[n])));
                                    //print out the estimated error. 
                                }
                                //Admittedly the meaning of the Big and Smol steps become quite convoluted when validating the
                                //adaptive methods. Thing of big as the higher order result, and smol as the lower order in this case. 

                                //Note: this will not produce an integer, but with proper data it will be close to an integer
                                //and the validation would be performed by rounding. 
                                //However one can also get errors if the results are too exact, roundoff error can ruin the calcluation. 
                                //Using larger step sizes usually removes that. 
                                for (int n = 0; n < numberOfEquations; n++) {
                                    yBigStep[n] = yBigStepHalf[n];
                                    ySmolSteps[n] = ySmolStepsHalf[n];
                                    //now undo what we just did so the corret values can be read.
                                }
                            }
                        } 

                    }
                }
                //Now that the step and double step have been taken, time to calculate some errors and see if we move 
                //on to the next step. 
                //First, from our parameters declared at the beginning, determine what our error limit is. 
                //Using GSL's version we frist estimate our errro based on what we know.
                if (i != 0) {
                    for (int n = 0; n<numberOfEquations; n++) {
                        errorEstimate[n] = sqrt((yBigStep[n] - ySmolSteps[n])*(yBigStep[n] - ySmolSteps[n]))* errorSafety;
                        //The 4/15 is taken from GSL's solver, a 'saftey factor' with unknown reasoning. 
                    }
                    errorEstimate[1+numberOfEquations] = sqrt((cpBigStep.rho - cpSmolSteps.rho)*(cpBigStep.rho - cpSmolSteps.rho))*errorSafety;
                    //find a way to generalize to any number of constants in the struct. 

                    double errorLimiter[numberOfEquations];
                    //since the definition of the error limiter uses a derivative, we cannot use it to limit the constant's error. 
                    //We originally had the error limiter set its own values. 
                    //GSL's formatting requries us to change this. 
                    system.function(currentPosition+step,ySmolSteps, errorLimiter, &cp);
                    //Now SmolSteps is used to set the errorLimiter. 
                    for (int n = 0; n<numberOfEquations; n++) {
                        errorLimiter[n] = absoluteErrorLimit + relativeErrorLimit*(ayErrorScaler*sqrt(ySmolSteps[n]*ySmolSteps[n]) + adyErrorScaler*step*sqrt(errorLimiter[n]*errorLimiter[n]));
                    }
                    
                    //The error limiter is set for every equation. Now we need to perform checks.

                    double ratioED = 0.0;
                    for (int n = 0; n<numberOfEquations; n++) { 
                        if (ratioED < errorEstimate[n]/errorLimiter[n]) {
                            ratioED = errorEstimate[n]/errorLimiter[n];
                            //pick out the largest of these ratios for use, every time. 
                        }
                    }

                    underError = false;
                    overError = false;
                    //make sure to set our values to false every loop. 

                    //these will be set to true when the condition is tripped. 
                    if (ratioED >  errorUpperTolerance) {
                        //If we are 10% (or whatever value is specified) over what the error we want is, adjust. 
                        overError = true;
                    } else if (ratioED <= errorLowerTolerance) {
                        //If we are 50% (or whatever value is specified) under what the error we want is, adjust. 
                        underError = true;
                    }

                    //If we have no trouble...
                    if (underError == false && overError == false) {
                        errorSatisfactory = true;
                    }
                    //...say that we're cleared to move to the next step. However, if one of them was triggered, we need to adjust. 
                    //in these cases we change the actual step size. 
                    //It is theoretically possible for both to be triggered on different equations. In that case, overError
                    //takes prescedent. We would rather have more accuracy than less in odd situations like that. 

                    //These if statements perform step adjustment if needed. Based on GSL's algorithm. 
                    else if (overError == true) {
                        step = step * scaleFactor * pow(ratioED,-1.0/butcher[rows-1-methodType][0]);
                        //printf("LOWER,%i %15.14e %15.14e %15.14e\n",i,currentPosition, step, ratioED);
                    } else { //if underError is true and overError is false is the only way to get here. The true-true situation is skipped.
                        step = step * scaleFactor * pow(ratioED,-1.0/(butcher[rows-1-methodType][0]+1));
                        errorSatisfactory = true;
                        //printf("UPPER,%i %15.14e %15.14e %15.14e\n",i,currentPosition, step, ratioED);
                    }

                    //Check to see if we're adjusting the step too much at once. 
                    //If we are, declare that we're done. 
                    if (step > maxStepAdjustment * originalStep) {
                        step = maxStepAdjustment * originalStep;
                        errorSatisfactory = true;
                    } else if (step < minStepAdjustment * originalStep){
                        step = minStepAdjustment * originalStep;
                        errorSatisfactory = true;
                    }

                    //We also declare some minium and maximum step conditions. 
                    if (step > absoluteMaxStep) {
                        step = absoluteMaxStep;
                        errorSatisfactory = true;
                    } else if (step < absoluteMinStep){
                        step = absoluteMinStep;
                        errorSatisfactory = true;
                    }

                    //With that, the step size has been changed. If errorSatisfactory should be false, it goes back and performs everything again
                    //with the new step size. 
                } else {
                    errorSatisfactory = true;
                    //We always want the *first* step to go through without change, often the first step is chosen for a specific reason. 
                    //In our work this generally came from a need to plot data sets against each other. 
                }
            }
            
            //Finally, we actually update the real answer. 
            for (int n = 0; n<numberOfEquations; n++) {
                if (methodType == 1) {
                    y[n]=yBigStep[n];
                } else {
                    y[n]=ySmolSteps[n];
                }
                //This check is required due to the way the butcher tables are stored. 
                //Probably a more efficient way to do this. 
            }

            if (underError == true) {
                currentPosition = currentPosition + originalStep;
                //if we had an underError and increased the step size, well, we kept the older points so we use that to update our current location.
            } else {
                currentPosition = currentPosition + step;
                //in any other case we use the new step. Even the case where the step wasn't actually changed. 
            }
            //After each step is calculated, print results. 
            //However, prior to printing we need to run our exception and constant evaluators one more time.
            exceptionHandler(currentPosition,y,&cp);
            constEval(currentPosition,y,&cp);
            //Since we've usually been running them on yInsert, the actual y and c values have not generally seen 
            //the restrictions applied here. 

            //Uncomment for live updates.
            //printf("Position:,\t%15.14e,\t",currentPosition);
            fprintf(fp, "Position:,\t%15.14e,\t",currentPosition);
            for (int n = 0; n < numberOfEquations; n++) {
                //printf("Equation %i:,\t%15.14e,\t",n, y[n]);
                fprintf(fp, "Equation %i:,\t%15.14e,\t",n, y[n]);
            }
            assignConstants(c,&cp); 
            for (int n = 0; n < numberOfConstants; n++) {
                //printf("Constant %i:,\t%15.14e,\t",n, c[n]);
                fprintf(fp, "Constant %i:,\t%15.14e,\t",n, c[n]);
            }
            //printf("\n");
            fprintf(fp,"\n");

            if (reportErrorEstimates == true) {
                //Print the error estimates we already have. 
                fprintf(fp, "Error Estimates:,\t");
                for (int n = 0; n < numberOfEquations; n++) {
                    fprintf(fp, "Equation %i:,\t%15.14e,\t",n,errorEstimate[n]);
                }
                for (int n = 0; n < numberOfConstants; n++) {
                    fprintf(fp, "Constant %i:,\t%15.14e,\t",n,errorEstimate[n+numberOfEquations]);
                }   
                fprintf(fp,"\n");
            }
            
            if (reportErrorActual == true) {
                //Now if we have an actual error to compare against with, there's some more work to do. 
                double yTruth[numberOfEquations];
                double cTruth[numberOfEquations];
                struct constantParameters cpTruth; 

                cpTruth.dimension = numberOfConstants;
                cpTruth.rho = cp.rho;
                knownQEval(currentPosition,yTruth);
                constEval(currentPosition,yTruth,&cpTruth);
                assignConstants(c,&cp); 
                assignConstants(cTruth,&cpTruth); 
                fprintf(fp, "Errors:,\t");
                for (int n = 0; n < numberOfEquations; n++) {
                    fprintf(fp, "Equation %i:,\t%15.14e,\t",n, yTruth[n]-y[n]);
                    fprintf(fp, "Truth:,\t%15.14e,\t",yTruth[n]);
                }
                for (int n = 0; n < numberOfConstants; n++) {
                    fprintf(fp, "Constant %i Error:,\t%15.14e,\t",n, cTruth[n]-c[n]);
                    fprintf(fp, "Truth:,\t%15.14e,\t",cTruth[n]);
                } 
                //printf("\n");
                fprintf(fp,"\n");
            }

            //printf("rho %10.9e \n", cp.rho);

            //And the very last thing we do in the loop is ask if we terminate it. 
            if (doWeTerminate(currentPosition, y, &cp) == 1) {
                i = SIZE;
            }
        }
    } else {
        //This loop is for the Adams-Bashforth method, which is implemented entirely differnetly from all RK methods.
        //As such it needs an entirely different algorithm. 
        for (int i = 0; i < SIZE; i++){ 
            //i represents how many steps have been taken. 0 is the initial condition, that is, the variable `bound`.

            //AB methods need to store values equal to the order of the method. 

            double yValues[numberOfEquations][(int)butcher[rows-1][0]];
            if (i == 0) {
                //first time initialization.
                for (int n = 0; n< numberOfEquations; n++) {
                    yValues[n][0] = y[n];
                    for (int m = 1; m < (int)butcher[rows-1][0]; m++) {
                        yValues[n][m] = 0; //these values shouldn't be used, but zero them anyway. 
                    } 
                }
            }
            //one for each equation, and then each equation also gets to store past values. 

            //validation is not available for AB methods since it doesn't really get to the steps of that order immediately. 

            struct constantParameters cpValues; 

            cpValues.dimension = numberOfConstants;
            cpValues.rho = cp.rho;
            //Find a way to genrealize this assignment to any number of items in the struct. 

            //This is normally where we would calulate the K values, but they are entirely unecessary here.

            double yInsert[numberOfEquations];
            //We also need an array for the inserted y-values for each equation. 
            //can we just use yValues directly? 
            //Most applications actually have the different yInsert values be independent, so 
            //if we knew the form of the equation we could simplify the code.
            //However, we need to make sure to always fill everything in case we have a system
            //of the form y'=f(u,y) u'=g(u,y)

            double dyOut[numberOfEquations];
            //GSL demands that we use two separate arrays for y and y', so here's y'. 

            struct constantParameters cpInsert; 

            cpInsert.dimension = numberOfConstants;
            cpInsert.rho = cp.rho;
            //Create an array to hold the constants we want.
            //Find way to generalize. 

            double xInsert = bound + step*i;

            //Check for any limitations on our results. 
            exceptionHandler(xInsert,yInsert,&cpInsert);

            //first, determine which row to use. 
            int currentRow;
            if (i < (int)butcher[rows-1][0]-1) {
                currentRow = (int)butcher[rows-1][0]-1-i;
                //basically, keep track of how many steps we actually have on offer to use. 
            } else {
                currentRow = 0;
                //the highest order part of the method is used. 
            }


            for (int m = (int)butcher[rows-1][0]-1; m >= 0; m--) {
                //we actually need m=0 in this case, the "present" is evaluated. 
                xInsert = bound + step*(i-m);
                //the "current locaiton" depends on how far in the past we are.
                for (int j = 0; j < numberOfEquations ; j++) {
                    yInsert[j] = yValues[j][m];
                }
                //Grab the correct yValues for the proper time/location. 

                //Check for any limitations on our results. 
                exceptionHandler(xInsert,yInsert,&cpInsert);

                //Evaluate the constants. 
                constEval(xInsert, yInsert,&cpInsert);

                //Now we actually evaluate the differential equations.
                system.function(xInsert, yInsert, dyOut, &cpInsert);

                //With that evaluation, we can change the value of y for each equation. 
                for (int n = 0; n< numberOfEquations; n++) {
                    y[n] = y[n] + step*butcher[n][m]*dyOut[n];
                }
                //Keep in mind this is procedural, y isn't right until all values of m have been cycled through. 
            }

            //At the END of every loop, we "shift" the values in the array down one space, that is, into the "past"
            //Present values are 0, previous step is 1, step before that is 2, etc. 
            for (int n = 0; n < numberOfEquations; n++) {
                for (int m = (int)butcher[rows-1][0] - 1; m > 0; m--) {
                    //possible error source here: loop indexing error, check later. 
                    yValues[n][m] = yValues[n][m-1];
                    //note that we start at the last column, m, and move the adjacent column to it. 
                    //This pushes off the value at the largest m value, since it's far enough in the past we no longer care.
                }
                yValues[n][0] = y[n]; 
                //present values update to what we just calculated. 
                //We have now completed stepping. 
            }            

            //Note that we specifically set ySmol to the value, not anything else. 
            //This is because we wish to avoid abusing if statements and would like to do the following only once.
            exceptionHandler(bound+step*(i+1),y,&cp);
            constEval(bound+step*(i+1),y,&cp); 
            //Check for exceptions and evaluate constants at our new values. 

            //Uncomment for live updates.
            //printf("Position:,\t%15.14e,\t",currentPosition);
            fprintf(fp, "Position:,\t%15.14e,\t",bound + step*(i+1));
            for (int n = 0; n < numberOfEquations; n++) {
                //printf("Equation %i:,\t%15.14e,\t",n, y[n]);
                fprintf(fp, "Equation %i:,\t%15.14e,\t",n, y[n]);
            }
            assignConstants(c,&cp); 
            for (int n = 0; n < numberOfConstants; n++) {
                //printf("Constant %i:,\t%15.14e,\t",n, c[n]);
                fprintf(fp, "Constant %i:,\t%15.14e,\t",n, c[n]);
            }
            //printf("\n");
            fprintf(fp,"\n");

            //can't report error esitmates.
            
            /*if (reportErrorActual == true) {
                //Now if we have an actual error to compare against with, there's some more work to do. 
                double yTruth[numberOfEquations];
                double cTruth[numberOfEquations];
                struct constantParameters cpTruth; 

                cpTruth.dimension = numberOfConstants;
                cpTruth.rho = cp.rho;
                knownQEval(currentPosition,yTruth);
                constEval(currentPosition,yTruth,&cpTruth);
                assignConstants(c,&cp); 
                assignConstants(cTruth,&cpTruth); 
                fprintf(fp, "Errors:,\t");
                for (int n = 0; n < numberOfEquations; n++) {
                    fprintf(fp, "Equation %i:,\t%15.14e,\t",n, yTruth[n]-y[n]);
                    fprintf(fp, "Truth:,\t%15.14e,\t",yTruth[n]);
                }
                for (int n = 0; n < numberOfConstants; n++) {
                    fprintf(fp, "Constant %i Error:,\t%15.14e,\t",n, cTruth[n]-c[n]);
                    fprintf(fp, "Truth:,\t%15.14e,\t",cTruth[n]);
                } 
                //printf("\n");
                fprintf(fp,"\n");
            }*/

            //printf("rho %10.9e \n", cp.rho);

            //And the very last thing we do in the loop is ask if we terminate it. 
            if (doWeTerminate(currentPosition, y, &cp) == 1) {
                i = SIZE;
            }
        }
    }
    //SECTION III: Analysis
    //Minor post-processing goes here. 
    //Anything advanced will need to be done in a data analysis program. 
    //We like to use matplotlib.

    // basic reference: https://www.tutorialspoint.com/cprogramming/c_file_io.htm
    // used to be a file converter here, now there isn't, we just close the file. 
    fclose(fp);

    printf("ODE Solver \"Odie\" V7 Shutting Down...\n");
    return 0;
}

// - GM, master of dogs.