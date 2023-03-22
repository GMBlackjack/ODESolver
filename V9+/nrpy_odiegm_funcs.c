#include "nrpy_odiegm_proto.c"

//This file contains the actual definitions for the funcitons outlined in nrpy_odiegm_proto.c

nrpy_odiegm_step *
nrpy_odiegm_step_alloc (const nrpy_odiegm_step_type * T, size_t dim)
{
  //Allocate the step "object", set all values, even those that may not be used. 
  nrpy_odiegm_step *s = (nrpy_odiegm_step *) malloc (sizeof (nrpy_odiegm_step));
  s->type = T;
  s->methodType = 1;
  s->adamsBashforthOrder = 0;
  s->rows = T->rows;
  s->columns = T->columns;
  //these last two assignments might be unecessary, but it will be convenient if this number
  //can be acessed at both levels. 
  if (T->rows == T->columns) {
    s->methodType = 0; //aka, normal RK-type method. 
  }
  if (T->rows == 19) {
    s->methodType = 2; //AB method. 
    s->adamsBashforthOrder = 4; //default order chosen, if user wants control they will specify elsewhere.
    //after allocation is run.  
  }

  s->yValues = (double *) malloc ((double)19.0 * dim * sizeof (double));
  //This here is the array used to store past values.
  //Only used for AB methods, but it still needs to be dynamically allocated. 
  //Having an adamsbashforthorder of 0 doesn't throw any errors, which is conveinent.

  return s;
}

nrpy_odiegm_evolve *
nrpy_odiegm_evolve_alloc (size_t dim)
{
  //allocate the evolve "object" and set all values, even those that may not be used.
  nrpy_odiegm_evolve *e = (nrpy_odiegm_evolve *) malloc (sizeof (nrpy_odiegm_evolve));
  e->y0 = (double *) malloc (dim * sizeof (double));
  e->yerr = (double *) malloc (dim * sizeof (double));
  //fill these with 0 just in case someone tries to allocate something. 
  for (int n = 0; n < dim; n++) {
    e->y0[n] = 0.0;
    e->yerr[n] = 0.0;
  }
  
  e->count = 0;
  e->last_step = 0.0; //By default we don't use this value. 
  e->bound = 0.0; //This will need to be adjusted to handle other starting positions. 
  e->currentPosition = e->bound;
  e->noAdaptiveTimestep = false; //We assume adaptive by default. 
  e->validate = false; //this one will need to be manually set if a user wants internal validation. 
  return e;
}

nrpy_odiegm_control *
nrpy_odiegm_control_y_new (double eps_abs, double eps_rel)
{
  //allocate the control "object." Unusual wording of function is due to us needing
  //a GSL replacement. 
  nrpy_odiegm_control *c = (nrpy_odiegm_control *) malloc (sizeof (nrpy_odiegm_control));
  c->absLim = eps_abs;
  c->relLim = eps_rel;

    c->scaleFactor = 0.9;
    c->errorSafety = 4.0/15.0;
    c->ayErrorScaler = 1.0;
    c->adyErrorScaler = 1.0;
    c->maxStepAdjustment = 5.0;
    c->minStepAdjustment = 0.2;
    c->absoluteMaxStep = 0.1;
    c->absoluteMinStep = 1e-10;
    c->errorUpperTolerance = 1.1;
    c->errorLowerTolerance = 0.5;
    //These are all the default values, virtually all responsible for adaptive timestep and 
    //error estimation.

  return c;
}

nrpy_odiegm_driver * nrpy_odiegm_driver_alloc_y_new (const nrpy_odiegm_system * sys,
                               const nrpy_odiegm_step_type * T,
                               const double hstart,
                               const double epsabs, const double epsrel)
{
    //Initializes an ODE driver "object" which contains all the "objets" above, making a system
    //that is prepared to evaluate a system of differential equations. 

    nrpy_odiegm_driver *state;
    state = (nrpy_odiegm_driver *) calloc (1, sizeof (nrpy_odiegm_driver)); //valgrind doesn't like this line. 
    const size_t dim = sys->dimension; 
    state->sys = sys;
    state->s = nrpy_odiegm_step_alloc (T, dim);

    state->e = nrpy_odiegm_evolve_alloc (dim);
    state->h = hstart; //the step size. 

    state->c = nrpy_odiegm_control_y_new (epsabs, epsrel);

  //There were functions here in GSL that assigned the driver to the objects contained in the driver.
  //We will not be doing that insanity. 

  return state;
}

//Memory freeing methods. 
void nrpy_odiegm_control_free (nrpy_odiegm_control * c)
{
  free (c);
}
void nrpy_odiegm_evolve_free (nrpy_odiegm_evolve * e)
{
  //free (e->dydt_out);
  //free (e->dydt_in); //these two are for jacobians, not used. 
  free (e->yerr);
  free (e->y0);
  free (e);
}
void nrpy_odiegm_step_free (nrpy_odiegm_step * s)
{ 
  free (s->yValues);
  free (s);
}
void nrpy_odiegm_driver_free (nrpy_odiegm_driver * state)
{
  //In most cases, this method should be called alone, calling the others would be redundant. 
  if (state->c)
    nrpy_odiegm_control_free (state->c);

  if (state->e)
    nrpy_odiegm_evolve_free (state->e);

  if (state->s)
    nrpy_odiegm_step_free (state->s);

  free (state);
}

//SECTION 2: Methods

//The actual stepping functions. 

//The goal is for these functions to be completely agnostic to whatever the user is doing, 
//they should always work regardless of the form of the system passed, the method passed, and even
//if the user does something dumb it shouldn't crash. It will spit out nonsense in those cases, though. 

int nrpy_odiegm_evolve_apply (nrpy_odiegm_evolve * e, nrpy_odiegm_control * c,
                             nrpy_odiegm_step * s,
                             const nrpy_odiegm_system * dydt, double *t,
                             double t1, double *h, double y[]) {
    //This is the big one, the function that ACTUALLY performs the timestep adjustment. 

    //First off, check if we're at the desired edge or not. 
    if (*t + *h > t1) {
        *h = t1 - *t;
        //If we're going past an endpoint we want, reduce the step size. 
        //Otherwise continue as normal. 
        //no need to stop the adaptive time step! If we need to increase the size, we
        //still report the smaller value, so it'll go through. 
        e->last_step = 1.0; //this is generally not used but the user might want it or something. 
    }

    //Gotta read in several things... improves readability.
    //Don't need a million arrows everywhere. 
    int numberOfEquations = (int)(dydt->dimension);
    double currentPosition = *t;
    double step = *h; 

    bool validate = e->validate;
    unsigned long int i = e->count;
    double bound = e->bound;
    bool noAdaptiveTimestep = e->noAdaptiveTimestep;

    int methodType = s->methodType; 
    int rows = s->type->rows;
    int columns = s->type->columns;
    int adamsBashforthOrder = s->adamsBashforthOrder;

    double absoluteErrorLimit = c->absLim;
    double relativeErrorLimit = c->relLim;
    double scaleFactor = c->scaleFactor;
    double errorSafety = c->errorSafety;
    double ayErrorScaler = c->ayErrorScaler;
    double adyErrorScaler = c->adyErrorScaler;
    double maxStepAdjustment = c-> maxStepAdjustment;
    double minStepAdjustment = c->minStepAdjustment;
    double absoluteMaxStep = c->absoluteMaxStep;
    double absoluteMinStep = c->absoluteMinStep;
    double errorUpperTolerance = c->errorUpperTolerance;
    double errorLowerTolerance = c->errorLowerTolerance;

    double yValues[numberOfEquations][adamsBashforthOrder];

    int counter = 0; //This counter is reused time and time again for sifting through memory
    //Allow me to express my dislike of void pointers. 

    //The following section only runs if we're using an AB method, otherwise it jumps over. 
    if (adamsBashforthOrder != 0) {
        if (i == 0) {
            //first time initialization of the yValues array for AB methods. 
            for (int n = 0; n< numberOfEquations; n++) {
                yValues[n][0] = y[n];
                for (int m = 1; m < adamsBashforthOrder; m++) {
                    yValues[n][m] = 0; //these values shouldn't be used, but zero them anyway. 
                } 
            }
        } else {
            //load values from known yValues if not first step for AB method. 
            for (int n = 0; n< numberOfEquations; n++) {
                for (int m = 0; m < adamsBashforthOrder; m++) {
                    yValues[n][m] = *((double *)(*s).yValues+counter); //Gotta fill in an array... joy...
                    counter++;
                    //This has to be done this way due to the array being passed as a void pointer. 
                } 
            }
        }
    }

    const nrpy_odiegm_step_type * stepType;
    stepType = s->type;

    counter = 0;
    if (methodType == 2) {
        rows = adamsBashforthOrder;
        columns = adamsBashforthOrder;
    }
    double butcher[rows][columns];
    //This is the butcher table that actually defines the method we use. 
    if (methodType != 2) { //If we aren't using AB method, just fill it without anything special. 
        for (int k=0; k < rows; k++) {
            for (int j = 0; j < columns; j++) {
                butcher[k][j] = *((double *)(*stepType).butcher+counter);
                counter++;
            }
        }
    } else { //If we ARE using an AB method, we need to construct it a little more carefully. 
        counter = counter + 19*(19-adamsBashforthOrder);
        //Every row has 19 elements, and we need to clear 19-order rows, leaving only the order behind. 
        for (int i=0; i < adamsBashforthOrder; i++) {
            counter = counter + 19-adamsBashforthOrder; 
            //for every row, clear the unneeded zeroes. 
            for (int j = 0; j < adamsBashforthOrder; j++) {
                butcher[i][j] = *((double *)(*stepType).butcher+counter);
                //This slowly counts through the array via complciated void pointer nonsense. 
                counter++;
            }
        }
        if (adamsBashforthOrder == 19) {
            butcher[adamsBashforthOrder-1][0] = 0.0;
            //Implementation artifact--we stored the order in the bottom left corner. 
            //but now we have to get rid of it so the algorithm doesn't freak out later.             
        }
    }

    if (methodType != 2) {
        //To use adaptive time-step, we need to store data at different step values:
        double yBigStep[numberOfEquations];
        double ySmolSteps[numberOfEquations];

        double yBigStepHalf[numberOfEquations];
        double ySmolStepsHalf[numberOfEquations];
        //These are here for validation only. We would rather not declare them at all, 
        //but they have to be declared outside the loop
        //so the information within can be stored through the iterations. 

        //One could argue that since the small steps will become our result 
        //we shouldn't declare it, however we are actually
        //NOT going to assign them to the actual answer y until we compare and run the adaptive
        //time-step algorithm. We might throw out all the data and need to run it again! 
        double errorEstimate[numberOfEquations];
        //even if we aren't limiting the constants, we can still report their error. 
        
        double originalStep = step;
        //We need to be able to refer to the original step so we can 
        //see if we're adjusting it too much at once. 

        //We rather explicitly do not actually take any steps until we confirm the error is below what we want.
        bool errorSatisfactory = false;
        bool underError = false;
        bool overError = false;
        //It's important to declare these outside the errorSatisfactory loop since to update the stepper we need to know
        //exactly what kind of step change we just did. 

        //This is a slapped together solution for indexing. 
        //Uses multiplication by 1 or 0 instead of an if statement on a bool. 
        //Should be more efficient. 
        int quickPatch = 1;
        if (methodType == 2) {
            quickPatch = 0;
        }
        //this constant removes certain components from consideraiton. 

        while (errorSatisfactory == false) {
            
            //All of the bellow values start off thinking they are the values from the 
            //previous step or initial conditions. 
            //We must reset them every time we return here.  
            for (int n = 0; n < numberOfEquations; n++) {
                yBigStep[n] = y[n];
                ySmolSteps[n] = y[n];
            } 
            for (int iteration = 1; iteration < 4; iteration++) {
                //So, we want to use Adaptive Timestep methodology. 
                //This will involve evaluating each step three times, 
                //In order to compare the evolution of two different step sizes and get an error estimate. 
                //Iteration 1 performs a normal step. 
                //Iteration 2 perofrms a half step.
                //Iteration 3 performs another half step after the previous one. 
                //Naturally the half-step results are reported as truth, 
                //but we get an error estimate from the difference
                //between the two values. 

                //For adaptive methods we only go through iteration 1 and 2
                
                //For AB method we only go through once, but do so with some additional operations. 

                if (i == 0 && iteration == 1 && validate == false && methodType == 0 && adamsBashforthOrder == 0) {
                    //don't take unecessary steps, if we are on the first step 
                    //and have no need for the large step, ignore it.
                    //Since we always want the first step to go through 
                    //don't bother calculating things we don't need. 
                    iteration = 2;
                    //This doesn't actually apply to inherently adaptive methods 
                    //since we cheat and do it in one iteration. 
                }

                double scale = 1.0;
                //this is the number we use to scale. It's either 1 or 1/2, 
                //depending on what size step we want. 
                int shift = 0;
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

                double K[rows-methodType*quickPatch][numberOfEquations];
                //These are the K-values that are required to evaluate RK-like methods. 
                //They will be determined based on the provided butcher table.
                //This is a 2D matrix since each diffyQ has its own set of K-values. 
                //Note that we subtract the method type from the row: adaptive RK butcher tables are larger. 

                //Since we'll be calling K while it's empty, even though there should be no errors due
                //to the way it's set up, let's go ahead and fill it with zeroes.
                for (int j = 0; j<rows-methodType*quickPatch; j++) {
                    for (int n = 0; n<numberOfEquations; n++) {
                        K[j][n]=0.0;
                    }
                } 

                double yInsert[numberOfEquations];
                //We also need an array for the inserted y-values for each equation. 

                double dyOut[numberOfEquations];
                //GSL demands that we use two separate arrays for y and y', so here's y'. 

                for (int j = 1; j < rows-methodType*quickPatch; j++) {
                    //Due to the way the butcher table is formatted, 
                    //start our index at 1 and stop at the end. 
                    double xInsert = currentPosition+shift*step*scale + butcher[j-1][0]*step*scale;
                    //xInsert does not change much for different tables, 
                    //just adjust the "step correction" term.
                    //xInsert is the same for every equation too.

                    for (int n = 0; n < numberOfEquations; n++) {
                        yInsert[n] = ySmolSteps[n];
                    } 
                    //Note that we are setting our buffer value, yInsert, to ySmolSteps. 
                    //This is because ySmolSteps is y at first, but we will need to evolve it forward
                    //two steps, so on the second small step this will be different. 
                    //(If using a method that requires that step, otherwise this is just a formality.)

                    for (int n = 1; n < columns; n++) {
                        //Once again, start at index of 1 rather than 0.
                        for (int q = 0; q < numberOfEquations; q++) {
                            yInsert[q] = yInsert[q] + butcher[j-1][n]*K[n][q];
                        }
                        //Each individual yInsert portion is dependent on one of the K values.
                        //K values are initially set to zero even though technically whenever 
                        //we would use an undeclared K-value the butcher table would have a zero.
                        //You know, just in case something goes wrong. 
                    }

                    //Now we actually evaluate the differential equations.
                    dydt->function(xInsert, yInsert, dyOut, dydt->params);
                    //yInsert goes in, dyOut comes out. Originally yInsert was overridden. 
                    //This no longer happens. 

                    for (int n = 0; n < numberOfEquations; n++) {
                        K[j][n] = step*scale*dyOut[n];
                        //Fill in the K-values we just calculated. 
                    } 
                }

                //Now that we have all the K-values set, we need to find 
                //the actual result in one final loop.
                for (int n = 0; n< numberOfEquations; n++) {
                    K[0][n] = ySmolSteps[n]; //The 0th spot in the K-values is reserved for holding the 
                    //final value while it's being calculated. 
                    for (int j = 1; j < columns; j++) {
                        K[0][n] = K[0][n] + butcher[rows-1-methodType*quickPatch][j]*K[j][n]; 
                        //This is where the actual approximation is finally performed. 
                    }
                    ySmolSteps[n] = K[0][n]; //Set ySmol to the new estimated value. 
                }
                //Note that we specifically set ySmol to the value, not anything else. 
                //This is because we wish to avoid abusing if statements.

                if (iteration == 1 || (i == 0 && validate == true && methodType == 1)) {
                    for (int n = 0; n<numberOfEquations; n++) {
                        yBigStep[n] = ySmolSteps[n];
                        ySmolSteps[n] = y[n];
                        //we still need to reset the value for SmolSteps on the first iteration
                        //no matter the type of method. 
                    }
                }
                //This only runs on the first iteration (or during validation), 
                //setting the big step to the right value
                //and resetting the small steps for when we actually use it. 
                //This odd structure exists purely for efficiency. 
                
                //If we are in an adaptive method situation, use that method and exit the iterations loop.
                if (methodType == 1) {
                    for (int n = 0; n< numberOfEquations; n++) {
                        K[0][n] = ySmolSteps[n]; //The 0th spot in the K-values is reserved 
                        //for holding the final value while it's being calculated. 
                        for (int j = 1; j < columns; j++) {
                            K[0][n] = K[0][n] + butcher[rows-1][j]*K[j][n]; 
                            //This is where the actual approximation is finally performed. 
                            //This time we use the bottom row, not the second to bottom row 
                            //(for adaptive methods)
                        }
                        ySmolSteps[n] = K[0][n]; //Set ySmol to the new estimated value. 
                    }
                    if (validate == true && i == 0 && iteration ==1) {
                        //do nothing fancy. 
                    } else {
                        iteration = 4; //break out after we get to the end, 
                        //we don't need to go any further. 
                        //We avoid this on the first iteration only if we want to validate. 
                    }
                }

                if (validate == true && i == 0) {
                    if (methodType == 0 && iteration == 2) {
                        //The normal validation for the normal RK methods. 

                        //Now that we've performed the approximation's first step at half the size, we can estimate the order.

                        //Create 2 arrays to hold the true values. 
                        double truthValidateBig[numberOfEquations];
                        double truthValidateSmol[numberOfEquations];
                        //Fill it with the true values.
                        dydt->trueFunction(bound+step,truthValidateBig);
                        dydt->trueFunction(bound+step*0.5,truthValidateSmol);
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
                    } else if (methodType == 1) {
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
                            dydt->trueFunction(bound+step,truthValidateBig);
                            dydt->trueFunction(bound+step*0.5,truthValidateSmol);
                            
                            //Then from this calculate the estimated errors.
                            for (int n = 0; n < numberOfEquations; n++) {
                                truthValidateBig[n] = (truthValidateBig[n] - yBigStepHalf[n]);
                                truthValidateSmol[n] = (truthValidateSmol[n] - yBigStep[n]);
                                //Now the validation steps contain their own errors, we can compare them.
                                printf("Order of ErrorA: %i\t%f\n",n, log2(truthValidateBig[n]/truthValidateSmol[n]));
                                //print out the estimated error. 
                            }
                            dydt->trueFunction(bound+step,truthValidateBig);
                            dydt->trueFunction(bound+step*0.5,truthValidateSmol);

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

                if (adamsBashforthOrder != 0) {
                    iteration = 4;
                    //We only iterate once for AB. Thus, break out if we are AB. 
                    for (int n = 0; n < numberOfEquations; n++) {
                        ySmolSteps[n] = yBigStep[n];
                    }
                }
            }
            //Now that the step and double step have been taken 
            //(or whatever alternative method is being used),
            //time to calculate some errors and see if we move on to the next step. 
            //First, from our parameters declared at the beginning, determine what our error limit is. 
            //Using GSL's version we frist estimate our errro based on what we know.
            if (i != 0 && adamsBashforthOrder == 0) {
                //Literally none of this is used for the AB method. 
                for (int n = 0; n<numberOfEquations; n++) {
                    errorEstimate[n] = sqrt((yBigStep[n] - ySmolSteps[n])*(yBigStep[n] - ySmolSteps[n]))* errorSafety;
                    //The 4/15 for errorSafety is taken from GSL's solver, a 'saftey factor' 
                    //with unknown reasoning. 
                }

                double errorLimiter[numberOfEquations];
                //since the definition of the error limiter uses a derivative, 
                //we cannot use it to limit the constant's error. 
                //We originally had the error limiter set its own values. 
                //GSL's formatting requries us to change this. 
                dydt->function(currentPosition+step,ySmolSteps, errorLimiter, dydt->params);
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

                counter = 0;
                for (int n = 0; n< numberOfEquations; n++) {
                    *((double *)(*e).yerr+counter) = errorEstimate[n]; //Gotta fill in an array... joy...
                    counter++;
                }
                //Notably this yerr is not used naturally by the program, it's set so the user can
                //read it if so desired. 

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
                if (noAdaptiveTimestep == false) {
                    //If we have no trouble...
                    if (underError == false && overError == false) {
                        errorSatisfactory = true;
                    }
                    //...say that we're cleared to move to the next step. 
                    //However, if one of them was triggered, we need to adjust. 
                    //In these cases we change the actual step size. 
                    //It is theoretically possible for both to be triggered on different equations. 
                    //In that case, overError takes prescedent. 
                    //We would rather have more accuracy than less in odd situations like that. 

                    //These if statements perform step adjustment if needed. Based on GSL's algorithm. 
                    else if (overError == true) {
                        step = step * scaleFactor * pow(ratioED,-1.0/butcher[rows-1-methodType*quickPatch][0]);
                    } else { //if underError is true and overError is false 
                        //is the only way to get here. The true-true situation is skipped.
                        step = step * scaleFactor * pow(ratioED,-1.0/(butcher[rows-1-methodType*quickPatch][0]+1));
                        errorSatisfactory = true;
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

                } else {
                    errorSatisfactory = true;
                    underError = false;
                }
                //With that, the step size has been changed. If errorSatisfactory is still false, 
                //it goes back and performs everything again with the new step size. 
            } else {
                errorSatisfactory = true;
                //We always want the *first* step to go through without change, 
                //often the first step is chosen for a specific reason. 
                //In our work this generally came from a need to plot data sets against each other. 
                //Also do this if we are using the AB method, as it has no error checks. 
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
            //There may be a more efficient way to do this. 
        }

        if (underError == true) {
            currentPosition = currentPosition + originalStep;
            //if we had an underError and increased the step size, 
            //well, we kept the older points so we use that to update our current location.
        } else {
            currentPosition = currentPosition + step;
            //in any other case we use the new step. Even the case where the step wasn't actually changed. 
        }

        //Before, the values were Printed here. This method no longer prints, 
        //printing is done outside any method. 

        if (adamsBashforthOrder > 0) {
            //At the END of every loop, we "shift" the values in the array "down" one space, 
            //that is, into the "past."
            //Present values are 0, previous step is 1, step before that is 2, etc. 
            for (int n = 0; n < numberOfEquations; n++) {
                for (int m = adamsBashforthOrder - 1; m > 0; m--) {
                    yValues[n][m] = yValues[n][m-1];
                    //note that we start at the last column, m, and move the adjacent column to it. 
                    //This pushes off the value at the largest m value, 
                    //since it's far enough in the past we no longer care.
                }
                yValues[n][0] = y[n]; 
                //present values update to what we just calculated. 
                //We have now completed stepping. 
            }  
        }
    } else {
        //This loop is for the Adams-Bashforth method, which is implemented 
        //entirely differnetly from all RK methods.
        //As such it needs an entirely different algorithm. 

        //This is normally where we would calulate the K values, but they are entirely unecessary here.

        double yInsert[numberOfEquations];
        //We also need an array for the inserted y-values for each equation. 

        double dyOut[numberOfEquations];
        //GSL demands that we use two separate arrays for y and y', so here's y'. 

        double xInsert; //This is generally going to be rather simple. 

        //first, determine which row to use in the AB butcher table. 
        int currentRow;
        if (i < adamsBashforthOrder-1) {
            currentRow = adamsBashforthOrder-1-i;
            //basically, keep track of how many steps we actually have on offer to use. 
        } else {
            currentRow = 0;
            //the highest order part of the method is used when we hit a certain step. 
        }

        for (int m = adamsBashforthOrder-currentRow-1; m >= 0; m--) {
            //we actually need m=0 in this case, the "present" is evaluated. 
            xInsert = bound + step*(i-m);
            //the "current locaiton" depends on how far in the past we are.
            for (int j = 0; j < numberOfEquations ; j++) {
                yInsert[j] = yValues[j][m];
            }
            //Grab the correct yValues for the proper time/location. 

            //Now we actually evaluate the differential equations.
            dydt->function(xInsert, yInsert, dyOut, dydt->params);

            //With that evaluation, we can change the value of y for each equation. 
            for (int n = 0; n< numberOfEquations; n++) {
                y[n] = y[n] + step*butcher[currentRow][m+currentRow]*dyOut[n];

            }
            //Keep in mind this is procedural, y isn't right until all values of m have been cycled through. 
        }

        //At the END of every loop, we "shift" the values in the array 
        //down one space, that is, into the "past"
        //Present values are 0, previous step is 1, step before that is 2, etc. 
        for (int n = 0; n < numberOfEquations; n++) {
            for (int m = adamsBashforthOrder-1; m > 0; m--) {
                yValues[n][m] = yValues[n][m-1];
                //note that we start at the last column, m, and move the adjacent column to it. 
                //This pushes off the value at the largest m value, 
                //since it's far enough in the past we no longer care.
            }
            yValues[n][0] = y[n]; 
            //present values update to what we just calculated. 
            //We have now completed stepping. 
        }         

        currentPosition = bound+step*(i+1);
        
        //AB Validation
        if (validate == true && i == adamsBashforthOrder && methodType == 2) {
            //validation is gonna be tricky for AB methods, as it requries exact values for everything.
            //Thus, let's create an array of exact values. 
            double yValidateValues[adamsBashforthOrder][numberOfEquations];
            //We will use this twice, once to do a sigle step, once to do a half step.
            //note the indeces are reversed from the usual yValues matrix. 
            
            for (int m = 0; m <adamsBashforthOrder; m++) {
                dydt->trueFunction(bound+step*(adamsBashforthOrder-m),yValidateValues[m]);
            }
            //This fills the entire yValidateValues array with exact values. 

            double yValidateBigStep[numberOfEquations];
            double yValidateBigStepResult[numberOfEquations];
            //We declare two variables for validation in regards to the larger step. 

            for (int n = 0; n< numberOfEquations; n++) {
                yValidateBigStepResult[n] = yValidateValues[0][n];
            }  
            //Take the actual true values for the present here. 

            for (int m = adamsBashforthOrder-currentRow-1; m >= 0; m--) {
                dydt->function(bound+step*(adamsBashforthOrder-m), yValidateValues[m], yValidateBigStep, dydt->params);
                //Evaluate the differential equations, store the result in yValidateBigStep
                for (int n = 0; n< numberOfEquations; n++) {
                    yValidateBigStepResult[n] = yValidateBigStepResult[n] + step*butcher[currentRow][m+currentRow]*yValidateBigStep[n];
                    //Perform the actual AB method. 
                }  
            }

            for (int m = 0; m <adamsBashforthOrder; m++) {
                dydt->trueFunction(bound+step*adamsBashforthOrder-step*0.5*m,yValidateValues[m]);
            }
            //Now set the exact values for half steps. The entire array is different
            //since reaching back in time by half steps alters positions of basically everything
            //except the present. 

            double yValidateSmolStep[numberOfEquations];
            double yValidateSmolStepResult[numberOfEquations];
            //Small values.

            for (int n = 0; n< numberOfEquations; n++) {
                yValidateSmolStepResult[n] = yValidateValues[0][n];
            }  
            //Insert the true values for the present. 

            //The below loop is the same as the previous one, except it makes sure to take half steps carefully.
            for (int m = adamsBashforthOrder-currentRow-1; m >= 0; m--) {
                dydt->function(bound+step*(adamsBashforthOrder)-step*0.5*m, yValidateValues[m], yValidateSmolStep, dydt->params);
                for (int n = 0; n< numberOfEquations; n++) {
                    yValidateSmolStepResult[n] = yValidateSmolStepResult[n] + step*0.5*butcher[currentRow][m+currentRow]*yValidateSmolStep[n];
                }   
            }
                
            for (int n = 0; n< numberOfEquations; n++) {
                dydt->trueFunction(bound+step*adamsBashforthOrder+step,yValidateBigStep);
                dydt->trueFunction(bound+step*adamsBashforthOrder+step*0.5,yValidateSmolStep);
            } 
            //Now that we have already used the "...Step" variables for calculation, they can now hold the truth values.       
            //Reduce, reuse, recycle after all!

            for (int n = 0; n < numberOfEquations; n++) {
                yValidateBigStep[n] = (yValidateBigStep[n] - yValidateBigStepResult[n]);
                yValidateSmolStep[n] = (yValidateSmolStep[n] - yValidateSmolStepResult[n]);
                //Now the validation steps contain their own errors, we can compare them.
                printf("Order of Error: %i\t%f\n",n, log2(yValidateBigStep[n]/yValidateSmolStep[n]));
                printf("And the other stuff...: %li %f %f\n",i, yValidateBigStep[n],yValidateSmolStep[n]);
                //print out the estimated error.
        }
    }
            
    }
    //Now we adjust any values that changed so everything outside the function can know it. 


    *h = step;
    *t = currentPosition;
    e->count = i+1;
    //Update things stored outside the function. 

    //update yValues, very important. We spent all that time shifting everything, we need to be able
    //to access it next time this function is called! 
    counter = 0;

    if (adamsBashforthOrder != 0) {
        //Put the new yValues back into the stored array. 
        for (int n = 0; n< numberOfEquations; n++) {
            for (int m = 0; m < adamsBashforthOrder; m++) {
                *((double *)(*s).yValues+counter) = yValues[n][m]; //Gotta fill in an array... joy...
                counter++;
            } 
        }
    }

    //in case the user needs it for some reason we also save the result to the evolve object.
    counter = 0;
    for (int n = 0; n< numberOfEquations; n++) {
        *((double *)(*e).y0+counter) = y[n]; //Gotta fill in an array... joy...
        counter++;
    }

    return 0;                      
}

int nrpy_odiegm_evolve_apply_fixed_step (nrpy_odiegm_evolve * e,
                                        nrpy_odiegm_control * con,
                                        nrpy_odiegm_step * step,
                                        const nrpy_odiegm_system * dydt,
                                        double *t, double h0,
                                        double y[]){
    //This method performs a single fixed time step. 
    e->noAdaptiveTimestep = true;
    nrpy_odiegm_evolve_apply(e, con, step, dydt, t, *t+h0, &h0, y);

    return 0;
}
int nrpy_odiegm_driver_apply (nrpy_odiegm_driver * d, double *t,
                             const double t1, double y[]){
    //Takes as many steps as requested at the driver level. 
    //Only really useful if you don't want to report anything until the end. Which. Sure.
    while (*t < t1) {
        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, d->sys, t, t1, &(d->h), y);
    }

    return 0;
}
int nrpy_odiegm_driver_apply_fixed_step (nrpy_odiegm_driver * d, double *t,
                                        const double h,
                                        const unsigned long int n,
                                        double y[]){
    //This just forces a fixed-step extrapolation. 
    d->e->noAdaptiveTimestep = true;
    nrpy_odiegm_driver_apply(d, t, h*(double)n, y);

    return 0;
}