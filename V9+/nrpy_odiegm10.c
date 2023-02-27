//#include "nrpy_odiegm.h"
#include "nrpy_odiegm_specific_methods.c"
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
/*struct constantParameters { 
    int dimension; //number that says how many we have. 
    double rho;
    //add more as necessary.  
};*/

int main()
{
    printf("Beginning ODE Solver \"Odie\" V10...\n");
    
    //CURRENTLY CHECKING THE HANDOVER.

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

    bool reportErrorEstimates = false;
    //enable to have every other line in "oData.txt" report the estimated erorr values.
    //AB methods cannot report error estimates. 
    //enable to have the actual errors reported. Produces junk data if there is no declared known function. 
    //BE WARNED: setting reporError (either kind) to true makes it print out all error data on another line, the file will have
    //to be read differently. 

    //ERROR PARAMETERS: Use these to set limits on the erorr. 
    double absoluteErrorLimit = 1e-10; //how big do we let the absolute error be?
    double relativeErrorLimit = 1e-10; //how big do we let the relative error be?
    //Note: the function for determining error creates a combined tester, it doesn't do either/or absolute or relative.
    //Constants should be set to 1 unless you know what you're dong, which I don't.

    //adaptive timestep constraints
    //We allow for some truly adaptable control to how the adaptive timestep works. All the paramaters can be adjusted. 

    bool noAdaptiveTimestep = false; //If you don't want to take adaptive timesteps, set this to true. 

    int adamsBashforthOrder = 13; //if using the AB method, specify which order you want.

    //We need to define a struct that can hold all possible constants. 
    struct constantParameters cp; 

    cp.dimension = numberOfConstants;
    //we'll set the actual parameters later. 

    nrpy_odiegm_system system = {diffyQEval,NULL,numberOfEquations,&cp};
    //This is the system of equations we solve.
    //The NULL is where the Jacobian would be, a holdover from GSL formatting, even though we never use the Jacobian. 

    //Now we set up the method. This is definitely a very roundabout way to do it,
    //And for now all the pointer nonsense is immediately undone, but here it is. 
    //All the butcher tables themselves are defined in butcher.c. 
    //We just need to create an object that gets them.
    const nrpy_odiegm_step_type * stepType;
    stepType = nrpy_odiegm_step_L6;
    //Here is where the method is actually set, by specific name since that's what GSL does. 

    const nrpy_odiegm_step_type * stepType2;
    stepType2 = nrpy_odiegm_step_L6;
    //this is a second step type "object" (struct) for hybridizing. 
    //Only used if the original type is AB.
    //Set to AB to use pure AB method. 

    nrpy_odiegm_driver *d;
    d = nrpy_odiegm_driver_alloc_y_new(&system, stepType,step, absoluteErrorLimit, relativeErrorLimit); //valgrind doesn't like this line.

    int methodType = 1;
    if (stepType->rows == stepType->columns) {
        methodType = 0; //aka, normal RK-type method. 
    } //no need for an else, we set it to 1 earlier to represent Adaptive methods. 
    //This integer is actually very useful as it can help us keep track of index changes between
    //the two types of methods! (We said, before we added another methodType that ruined that).

    //Technically the above sets the type. What we do now is basically undo that pointer nonsense. 
    int rows = stepType->rows;
    int columns = stepType->columns;
    //Since we know the dimension, we can now fill an actual butcher array. 
    if (rows == 19) {
        //We're using the Adams method, special table construction required. 
        methodType = 2;
        noAdaptiveTimestep = true;
        reportErrorEstimates = false;
        if (stepType2 == nrpy_odiegm_step_AB) {
            rows = adamsBashforthOrder;
            columns = adamsBashforthOrder;
        } else {
            rows = stepType2->rows;
            columns = stepType2->columns;
        }
    } else {
        adamsBashforthOrder == 0;
    }
    double butcher[rows][columns];
    double butcher2[adamsBashforthOrder][adamsBashforthOrder];
    //We have two tables. The second one is only used for hybrid AB methods, it remains empty otherwise. 
    //In pure AB method, both the tables will be identical, but in hybrid, one will be RK (butcher) and the other AB (butcher2)

    //How to get array size no longer needed, was part of the array reading process. 

    if (validate == true) {
        printf("Method Order: %i. \nOrder of Error should be near to or larger than Method Order + 1.\n",stepType->order);
        printf("If not, try a larger step size, roundoff error may be interfering.\n");
    } else {
        printf("Method Order: %i.\n",stepType->order);
    }
    //If validation is not needed, we don't care about the Order of the Error. 
    //Currently not working, order of method will need to be reported another way. For now, though, valgrind it. 

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
    //When we adjust the step size it is not possible to algorithmically determine our position. 
    //Thus this variable is required.

    //This here sets the initial conditions as declared in getInitialCondition()
    getInitialCondition(y); 

    //This evaluates any constants that might be needed for evaluating the actual differnetial equations. 
    constEval(currentPosition, y, &cp);

    FILE *fp2;
    fp2 = fopen("ooData.txt","w");
    //This is where we test the GSL methods without destroying anything else. 

    //we assume that there are no exceptions to check for on initial conditions. 

    //SECTION II: The Loop


    //Open the file we'll be writing data to. 
    //First though, let's print out our initial data. The print function needs to be adaptable to any size of data. 
    //We can do this with multiple print functions and just not adding the newline character until we're done.
    //We print both to console and to the file for the initial conditions, but later only print to file.
    //First, print the location we are at. 
    printf("INITIAL: Position:,\t%f,\t",bound);
    //Second, go through and print the result for every single equation in our system.
    for (int n = 0; n < numberOfEquations; n++) {
        printf("Equation %i:,\t%15.14e,\t",n, y[n]);
    }
    //Third, print out desired constants.
    assignConstants(c,&cp);    
    for (int n = 0; n < numberOfConstants; n++) {
        printf("Constant %i:,\t%15.14e,\t",n, c[n]);
    }
    //Lastly, the newline character. 
    printf("\n");
    //Comma delimiters are printed to the file so it can be converted to .csv with ease. 

    /*if (reportErrorEstimates == true) {
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
    }*/

        //Secondary Test:

    //END GSL TEST.

    //Keep in mind that if both are set to true, both lines are printed. This is to ensure the genreated file cycles through 
    //rows of output identically at the start and throughout the program. 

    //The below array is only used for AB methods. 
    //However it still has to be declared outside the loop since it's used in all sections when it is needed. 
    //It stores the values attained in the present and the past up to the order of the AB method. 
    //Notably when not using an AB method adamsBashforthOrder == 0 so this shouldn't take up space. 

    //This loop fills out all the data.
    //It takes a provided butcher table and executes the method stored within. Any table should work.  

    int holder; 

    for (int i = 0; i < SIZE; i++){
        //And here we need some work... 

        if (methodType == 2 && i < adamsBashforthOrder && stepType2 != nrpy_odiegm_step_AB) {
            d->s->type = stepType2;
            d->s->rows = stepType2->rows;
            d->s->columns = stepType2->columns;
            d->s->methodType = 0;
            d->s->adamsBashforthOrder = adamsBashforthOrder;
            d->e->noAdaptiveTimestep = true;
        } else if (stepType != stepType2 && methodType == 2 && i >= adamsBashforthOrder) {
            d->s->type = stepType;
            d->s->rows = stepType->rows;
            d->s->columns = stepType->columns;
            d->s->methodType = 2;
            d->s->adamsBashforthOrder = adamsBashforthOrder;
            d->e->noAdaptiveTimestep = true;
        }

        holder = nrpy_odiegm_evolve_apply(d->e, d->c, d->s, &system, NULL, 0.0, &step, y);


            //Uncomment for live updates.
            //printf("Position:,\t%15.14e,\t",currentPosition);
            fprintf(fp2, "Position:,\t%15.14e,\t",d->e->currentPosition);
            for (int n = 0; n < numberOfEquations; n++) {
                //printf("Equation %i:,\t%15.14e,\t",n, y[n]);
                fprintf(fp2, "Equation %i:,\t%15.14e,\t",n, y[n]);
            }

            constEval(currentPosition, y, &cp);
            assignConstants(c,&cp);  


            for (int n = 0; n < numberOfConstants; n++) {
                //printf("Constant %i:,\t%15.14e,\t",n, c[n]);
                fprintf(fp2, "Constant %i:,\t%15.14e,\t",n, c[n]);
            }
            //printf("\n");
            fprintf(fp2,"\n");

            

            if (doWeTerminate(d->e->currentPosition, y, &cp) == 1) {
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