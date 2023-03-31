#include "nrpy_odiegm_funcs.c" //nrpy_odiegm itself.
#include "nrpy_odiegm_user_methods.c" //user-dependent functions. 
//TODO:
//Jupyter notebook Adjusents. 

//This file is technically not part of Odie, it is just an example implementation.
//However, it is exceptionally versatile, and can be used to run virtually 
//every system of differential equations desired.
//However, since it's general, it trades off some efficiency to handle
//every type of method. If efficiency is desired, it is recommended that the user 
//make a custom main(). 

int main()
{
    printf("Beginning ODE Solver \"Odie\" V10...\n");

    //SECTION I: Preliminaries
    //Before the program actually starts, variables need to be created
    //and set, as well as the functions chosen. 
    //The system of differential equations can be found declared in diffyQEval()
    //in nrpy_odiegm_user_methods.c

    double step = 0.00001; //the "step" value. Initial step if using an adaptive method.
    double currentPosition = 0.0; //where the boundary/initial condition is. Same for every equation in the system.
    int numberOfEquations = 4; //How many equations are in our system?
    int numberOfConstants = 1; //How many constants do we wish to separately evaluate and report? 
    //If altering the two "numberOf" ints, be careful it doesn't go over the actual number 
    //and cause an overflow in the functions in user_methods
    const int SIZE = 100000; //How many steps are we going to take? 
    //This is the default termination condition. 
    int adamsBashforthOrder = 4; //if using the AB method, specify which order you want.
    //If we are not using the AB method this is set to 0 later automatically. 4 by default. 
    bool noAdaptiveTimestep = false; //Sometimes we just want to step forward uniformly 
    //without using GSL's awkward setup. False by default. 

    bool reportErrorActual = false;
    bool reportErrorEstimates = false;
    //AB methods do not report error estimates. 
    //BE WARNED: setting reporError (either kind) to true makes it print out all error data on another line,
    //the file will have to be read differently. 

    //ERROR PARAMETERS: Use these to set limits on the erorr. 
    double absoluteErrorLimit = 1e-14; //how big do we let the absolute error be?
    double relativeErrorLimit = 1e-14; //how big do we let the relative error be?
    //Note: there are a lot more error control numbers that can be set inside the control "object" d->c.

    char fileName[] = "ooData.txt"; //Where do you want the data to print?

    //Now we set up the method. 
    const nrpy_odiegm_step_type * stepType;
    stepType = nrpy_odiegm_step_RK4;
    //Here is where the method is actually set, by specific name since that's what GSL does. 

    const nrpy_odiegm_step_type * stepType2;
    stepType2 = nrpy_odiegm_step_RK4;
    //this is a second step type "object" (struct) for hybridizing. 
    //Only used if the original type is AB.
    //Set to AB to use pure AB method. 

    //AFTER THIS POINT THERE SHOULD BE NO NEED FOR USER INPUT, THE CODE SHOULD HANDLE ITSELF. 

    //We need to define a struct that can hold all possible constants. 
    struct constantParameters cp; 
    cp.dimension = numberOfConstants;
    //we'll set the actual parameters later. 
    //Do note that cp itself needs to be declared in constantParameters in user_methods manually.
    //The methods that work with it need to be declared as well. 

    nrpy_odiegm_system system = {diffyQEval,knownQEval,numberOfEquations,&cp};
    //This is the system of equations we solve.
    //The second slot was originally the Jacobian in GSL, but we use it to pass a 
    //validation equation that may or may not be used.

    nrpy_odiegm_driver *d;
    d = nrpy_odiegm_driver_alloc_y_new(&system, stepType, step, absoluteErrorLimit, relativeErrorLimit); 
    //This is the "object" that runs everything, contains every needed varaible, etc. 
    //Basically the master of the whole thing, hence why it's called the "driver"
    //Contains three major sub-objects besides the step type. 
    //c is the controller, which is primarily used to store adaptive timestep values. 
    //s is the step, which has the step type in it, but also parameters that describe the steps themselves.
    //e is the evolver, which actually performs the update when it is requested. 

    int methodType = 1;
    if (stepType->rows == stepType->columns) {
        methodType = 0; //aka, normal RK-type method. 
    } //no need for an else, we set it to 1 earlier to represent Adaptive methods. 
    //This integer is actually very useful as it can help us keep track of index changes between
    //the two types of methods! (We said, before we added another methodType that ruined that).

    //Technically the above sets the type. What we do next is basically undo the pointer nonsense. 
    //Since we know the dimension, we can now fill an actual butcher array. 
    if (stepType->rows == 19) { 
        methodType = 2;
    } else {
        adamsBashforthOrder = 0;
    }
    d->s->adamsBashforthOrder = adamsBashforthOrder;
    d->e->noAdaptiveTimestep = noAdaptiveTimestep;
    //based on what type of method we are using, we adjust some parameters within the driver.

        if (methodType == 2) {
            printf("Method Order: %i.\n",adamsBashforthOrder);
        } else {
            printf("Method Order: %i.\n",stepType->order);            
        }
    
    double y[numberOfEquations];
    //These variables temporarily store the values calculated before they are 
    //printed to the output file and forgotten.
    //y contains the values of the actual equations and their derivatives. 
    //Each array only holds values at one evaluation point, but one for each Equation.

    double c[numberOfConstants];
    //c is just used to hold any constants we wish to report. 
    //You'd think that, since we have the constants in a struct, we can avoid declaring this.
    //No. Not as far as we can tell, anyway. Structs are a pain to iterate through,
    //and we can't know what form the user is going to hand us the struct in. 

    //This here sets the initial conditions as declared in getInitialCondition()
    getInitialCondition(y); 
    constEval(currentPosition, y,&cp);

    FILE *fp2;
    fp2 = fopen(fileName,"w");
    printf("Printing to file '%s'.\n",fileName);

    //Open the file we'll be writing data to. 
    //Before continuing, let's print out our initial data. 
    //The print function needs to be adaptable to any size of data. 
    //We can do this with multiple print functions and just 
    //not adding the newline character until we're done.
    //We print both to console and to the file for the initial conditions, but later only print to file.
    //First, print the location we are at. 
    printf("INITIAL: Position:,\t%f,\t",currentPosition);
    fprintf(fp2, "Position:,\t%15.14e,\t",currentPosition);
    //Second, go through and print the result for every single equation in our system.
    for (int n = 0; n < numberOfEquations; n++) {
        printf("Equation %i:,\t%15.14e,\t",n, y[n]);
        fprintf(fp2, "Equation %i:,\t%15.14e,\t",n, y[n]);
    }
    //Third, print out desired constants.
    assignConstants(c,&cp);    
    for (int n = 0; n < numberOfConstants; n++) {
        printf("Constant %i:,\t%15.14e,\t",n, c[n]);
        fprintf(fp2, "Constant %i:,\t%15.14e,\t",n, c[n]);
    }
    //Lastly, the newline character. 
    printf("\n");
    fprintf(fp2,"\n");
    //Comma delimiters are printed to the file so it can be converted to .csv with ease. 

    if (reportErrorEstimates == true) {
        //In order to keep things neat and regular in the file, print a first line of errors. 
        //Even though by necessity all of them must be zero. 
        fprintf(fp2, "Errors Estimates:,\t");
        for (int n = 0; n < numberOfEquations; n++) {
            fprintf(fp2, "Equation %i:,\t0.0,\t",n);
        }
        for (int n = 0; n < numberOfConstants; n++) {
            fprintf(fp2, "Constant %i:,\t0.0,\t",n);
        }   
        fprintf(fp2,"\n");
    }
    
    if (reportErrorActual == true) {
        //In order to keep things neat and regular in the file, print a first line of errors. 
        //Even though by necessity all of them must be zero. 
        fprintf(fp2, "Errors:,\t");
        for (int n = 0; n < numberOfEquations; n++) {
            fprintf(fp2, "Equation %i:,\t0.0,\t",n);
            fprintf(fp2, "Truth:,\t%15.14e,\t",y[n]);
        }
        for (int n = 0; n < numberOfConstants; n++) {
            fprintf(fp2, "Constant %i:,\t0.0,\t",n);
            fprintf(fp2, "Truth:,\t%15.14e,\t",c[n]);
        }   
        //printf("\n");
        fprintf(fp2,"\n");
    }


    //SECTION II: The Loop

    //This loop fills out all the data.
    //It takes a provided butcher table and executes the method stored within. Any table should work.  

    for (int i = 0; i < SIZE; i++){
        
        //Hybrid Methods require some fancy footwork, hence the if statements below. 
        if (methodType == 2 && i == 0 && stepType2 != nrpy_odiegm_step_AB) {
            d->s->type = stepType2;
            d->s->rows = stepType2->rows;
            d->s->columns = stepType2->columns;
            d->s->methodType = 0;
            d->s->adamsBashforthOrder = adamsBashforthOrder;
            d->e->noAdaptiveTimestep = true;
        } else if (stepType != stepType2 && methodType == 2 && i == adamsBashforthOrder) {
            d->s->type = stepType;
            d->s->rows = stepType->rows;
            d->s->columns = stepType->columns;
            d->s->methodType = 2;
            d->s->adamsBashforthOrder = adamsBashforthOrder;
            d->e->noAdaptiveTimestep = true;
        }

        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, &system, &currentPosition, currentPosition+step, &step, y);
        
        exceptionHandler(currentPosition,y);
        constEval(currentPosition,y,&cp);
        assignConstants(c,&cp);
        //These lines are to make sure the constant updates. 
        //And exception constraints are applied.  

        //Uncomment for live updates.
        //printf("Position:,\t%15.14e,\t",currentPosition);
        fprintf(fp2, "Position:,\t%15.14e,\t",currentPosition);
        for (int n = 0; n < numberOfEquations; n++) {
            //printf("Equation %i:,\t%15.14e,\t",n, y[n]);
            fprintf(fp2, "Equation %i:,\t%15.14e,\t",n, y[n]);
        }

        for (int n = 0; n < numberOfConstants; n++) {
            //printf("Constant %i:,\t%15.14e,\t",n, c[n]);
            fprintf(fp2, "Constant %i:,\t%15.14e,\t",n, c[n]);
            //printf("Constant %i:,\t%15.14e %15.14e,\n",n, c[n], y[n]);
        }
        //printf("\n");
        fprintf(fp2,"\n");

        if (reportErrorEstimates == true) {
            //Print the error estimates we already have. 
            fprintf(fp2, "Error Estimates:,\t");
            for (int n = 0; n < numberOfEquations; n++) {
                fprintf(fp2, "Equation %i:,\t%15.14e,\t",n,*(d->e->yerr)); //find a way to grab the error. 
            }
            //constant estimates not reported, only diffyQ values. 
            fprintf(fp2,"\n");
        }
            
        if (reportErrorActual == true) {
            //Now if we have an actual error to compare against with, there's some more work to do. 
            double yTruth[numberOfEquations];
            double cTruth[numberOfConstants];
            struct constantParameters cpTruth; 
            //True values for everything we compare with.
            
            knownQEval(currentPosition,yTruth);
            constEval(currentPosition,yTruth,&cpTruth);

            assignConstants(c,&cp); 
            assignConstants(cTruth,&cpTruth);
 
            fprintf(fp2, "Errors:,\t");
            for (int n = 0; n < numberOfEquations; n++) {
                fprintf(fp2, "Equation %i:,\t%15.14e,\t",n, yTruth[n]-y[n]);
                fprintf(fp2, "Truth:,\t%15.14e,\t",yTruth[n]);
            }
            for (int n = 0; n < numberOfConstants; n++) {
                fprintf(fp2, "Constant %i Error:,\t%15.14e,\t",n, cTruth[n]-c[n]);
                fprintf(fp2, "Truth:,\t%15.14e,\t",cTruth[n]);
            } 
            //printf("\n");
            fprintf(fp2,"\n");
        }

        if (doWeTerminate(currentPosition, y, &cp) == 1) {
            i = SIZE;
        }
    }

    //SECTION III: Analysis
    //Minor post-processing goes here. 
    //Anything advanced will need to be done in a data analysis program. 
    //We like to use matplotlib for python.

    // basic reference: https://www.tutorialspoint.com/cprogramming/c_file_io.htm
    // used to be a file converter here, now there isn't, we just close the file. 

    fclose(fp2);

    nrpy_odiegm_driver_free(d);
    //MEMORY SHENANIGANS

    printf("ODE Solver \"Odie\" V10 Shutting Down...\n");
    return 0;
}

// - GM, master of dogs.