#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

struct constantParameters { 
    int dimension; //number that says how many we have. 
    double rho;
    //add more as necessary.  
};

typedef struct {
    int (*function) (double x, const double y[], double dydx[], void *params);
    //The function passed to this struct contains the definitions of the differnetial equations. 
    int (*jacobian) (double t, const double y[], double *dfdy, double dfdt[],
                   void *params); //The jacobian is a holdover from GSL, it will not be used in this program,
                   //Always pass NULL to this part. 
                   //If we were not making a drop-in replacement, we would just have a three-paraemeter struct.
    size_t dimension; //For storing how big our system of equations is. 
    //Just pass it an int, usually. 
    void *params; //For storing extra constants needed to evaluate the functions. 
    //params->dimension stores how many there are. 
} nrpy_odiegm_system;

typedef struct {
    //unlike with the system object above, this step_type object does not need
    //to match GSL's form explicitly, it just needs to define the method.
    //Since we are using RK-type methods exclusively, this is just a
    int rows; 
    int columns; //since we are passing a void pointer to do this, we need a way
    //to know how large it is in the end.
    //void *(*alloc) (size_t dim); //Adding this breaks the array, we can't have it. 
    void *butcher;
    //Make sure to put this at the end of the struct though,
    //in case we add more parts to it, nonspecific arrays must be the last element.

    //two of these step_type objects might be needd at once, depending on implementation. 
} nrpy_odiegm_step_type;

typedef struct {
  const nrpy_odiegm_step_type *type; 
  int methodType; //can be calculated from above, but is easier to just store. 0,1,2 values. 
  int adamsBashforthOrder; //order if an AB method is used.
  size_t dimension;
  void *state;
} nrpy_odiegm_step;

typedef struct {
    double absLim;
    double relLim;
    double errorSafety;
    double ayErrorScalar;
    double adyErrorScalar;
    double maxStepAdjustment;
    double minStepAdjustment;
    double absoluteMaxStep;
    double absoluteMinStep;
    //We added these ourselves. 

    //const nrpy_odiegm_control_type *type;
    void *state;
  
} nrpy_odiegm_control;

typedef struct
{
  size_t dimension;
  double *y0; //the y array that stores the values. 
  double *yerr;
  double *dydt_in;
  double *dydt_out;
  double last_step;
  double bound; //the point at which we started is sometimes important. 
  unsigned long int count;
  unsigned long int failed_steps;
  bool validate;
  //const nrpy_odiegm_driver *driver;
  //NO SELF REFERENCING LOOPS! BAD! 
  void *yValues; //The extremely funky parameter that hides a 2D array, used when
  //the past steps are important for AB method. 
} nrpy_odiegm_evolve;

typedef struct {
    const nrpy_odiegm_system *sys; /* ODE system */
    nrpy_odiegm_step *s;           /* stepper object */
    nrpy_odiegm_control *c;        /* control object */
    nrpy_odiegm_evolve *e;         /* evolve object */
    double h;                     /* step size */
    double hmin;                  /* minimum step size allowed */
    double hmax;                  /* maximum step size allowed */
    unsigned long int n;          /* number of steps taken; i */
    unsigned long int nmax;       /* Maximum number of steps allowed; SIZE*/
} nrpy_odiegm_driver;

//Allocation Functions below. 

nrpy_odiegm_step *
nrpy_odiegm_step_alloc (const nrpy_odiegm_step_type * T, size_t dim)
{
  nrpy_odiegm_step *s = (nrpy_odiegm_step *) malloc (sizeof (nrpy_odiegm_step));
  s->type = T;
  s->dimension = dim;
  //s->state = s->type->alloc (dim);
  //We will see if we can avoid setting the state. 
  return s;
}

nrpy_odiegm_evolve *
nrpy_odiegm_evolve_alloc (size_t dim)
{
  nrpy_odiegm_evolve *e = (nrpy_odiegm_evolve *) malloc (sizeof (nrpy_odiegm_evolve));
  e->y0 = (double *) malloc (dim * sizeof (double));
  e->yerr = (double *) malloc (dim * sizeof (double));
  e->dydt_in = (double *) malloc (dim * sizeof (double));
  e->dydt_out = (double *) malloc (dim * sizeof (double));
  e->dimension = dim;
  e->count = 0;
  e->failed_steps = 0;
  e->last_step = 0.0;
  //e->driver = NULL; 
  //LOOP! BAD! *newspaper bap to snoot*
  return e;
}

nrpy_odiegm_control *
nrpy_odiegm_control_y_new (double eps_abs, double eps_rel)
{
  //return gsl_odeiv2_control_standard_new (eps_abs, eps_rel, 1.0, 0.0);
  //This was all that was in the original function. I once again wonder what GSL was thinking.
  nrpy_odiegm_control *c = (nrpy_odiegm_control *) malloc (sizeof (nrpy_odiegm_control));
  c->absLim = eps_abs;
  c->relLim = eps_rel;
  return c;
}

nrpy_odiegm_driver * nrpy_odiegm_driver_alloc_y_new (const nrpy_odiegm_system * sys,
                               const nrpy_odiegm_step_type * T,
                               const double hstart,
                               const double epsabs, const double epsrel)
{
    //Initializes an ODE driver system with control object of type y_new.

    nrpy_odiegm_driver *state;
    state = (nrpy_odiegm_driver *) calloc (1, sizeof (nrpy_odiegm_driver));
    const size_t dim = sys->dimension;
    state->sys = sys;
    state->s = nrpy_odiegm_step_alloc (T, dim);
    state->e = nrpy_odiegm_evolve_alloc (dim);
    state->h = hstart;
    state->hmin = 0.0;
    state->hmax = 1.0;
    state->nmax = 0;
    state->n = 0;
    state->c = NULL;
  if (epsabs >= 0.0 && epsrel >= 0.0)
    {
      state->c = nrpy_odiegm_control_y_new (epsabs, epsrel);
      //Otherwise we're not using a method that needs step size control. 
    }
  //There were functions here that assigned the driver to the objects contained in the driver.
  //We will not be doing that nonsense. 

  return state;
}

//Memory freeing methods. 
//Only necessary since we want a drop in replacement.
void nrpy_odiegm_control_free (nrpy_odiegm_control * c)
{
  //RETURN_IF_NULL (c);
  //c->type->free (c->state);
  free (c);
}
void nrpy_odiegm_evolve_free (nrpy_odiegm_evolve * e)
{
  //RETURN_IF_NULL (e);
  free (e->dydt_out);
  free (e->dydt_in);
  free (e->yerr);
  free (e->y0);
  free (e);
}
void nrpy_odiegm_step_free (nrpy_odiegm_step * s)
{
  //RETURN_IF_NULL (s);
  //s->type->free (s->state);
  //free (s->type);
  //These look to be unneeded. 
  free (s);
}
void nrpy_odiegm_driver_free (nrpy_odiegm_driver * state)
{
  if (state->c)
    nrpy_odiegm_control_free (state->c);

  if (state->e)
    nrpy_odiegm_evolve_free (state->e);

  if (state->s)
    nrpy_odiegm_step_free (state->s);

  free (state);
}

//The actual stepping functions. 
int nrpy_odiegm_evolve_apply (nrpy_odiegm_evolve * e, nrpy_odiegm_control * con,
                             nrpy_odiegm_step * step,
                             const nrpy_odiegm_system * dydt, double *t,
                             double t1, double *h, double y[]) {
    //This is the big one, the function that ACTUALLY performs the timestep adjustment. 
                //If we're not doing Adams-Bashforth, do the "normal" loop for RK-like methods.
                //For now, we just want to get THIS part working, AB is not our concern, but we will keep its checks in for now. 
                //OR do the "normal" loop if using a hybrid AB approach, until we hit the AB order.
                //i represents how many steps have been taken. 0 is the initial condition, that is, the variable `bound`.

                //Gotta read in a few things...

                double numberOfConstants = (*(struct constantParameters*)dydt->params)).dimension
                struct constantParameters cp;
                cp.dimension = numberOfConstants;


                //To use adaptive time-step, we need to store data at different step values:
                double yBigStep[dydt->dimension];
                double ySmolSteps[dydt->dimension];

                double yBigStepHalf[dydt->dimension];
                double ySmolStepsHalf[dydt->dimension];
                //These are here for validation only. We would rather not declare them at all, but they have to be declared outside the loop
                //So the information within can be stored through the iterations. 

                struct constantParameters cpBigStep; 

                updateConstantParameters(&cpBigStep, &cp);

                struct constantParameters cpSmolSteps; 

                updateConstantParameters(&cpSmolSteps, &cp);
                //This is the rather ugly way to assign the values of the constants we care about. 
                //Notably this is NOT GENERALIZED, and it NEEDS TO BE.

                //One could argue that since the small steps will become our result we shouldn't declare it, however we are actually
                //NOT going to assign them to the actual answer y until we compare and run the adaptive
                //time-step algorithm. We might throw out all the data and need to run it again! 
                double errorEstimate[dydt->dimension+numberOfConstants];
                //even if we aren't limiting the constants, we can still report their error. 
                
                double originalStep = step;
                //We need to be able to refer to the original step so we can see if we're adjusting it too much at once. 

                //We rather explicitly do not actually take any steps until we confirm the error is below what we want.
                bool errorSatisfactory = false;
                bool underError = false;
                bool overError = false;
                //It's important to declare these outside the errorSatisfactory loop since to update the stepper we need to know
                //exactly what kind of step change we just did. 

                    //This is a slapped together solution for indexing. 
                    //Find a better way, perhaps? 
                int quickPatch = 1;
                if (methodType == 2) {
                    quickPatch = 0;
                }
                //this constant just removes certain components from consideraiton. 

                while (errorSatisfactory == false) {
                    
                    //All of the bellow values start off thinking they are the values from the previous step or initial conditions. 
                    //We must reset them every time we return here.  
                    for (int n = 0; n < dydt->dimension; n++) {
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
                        
                        //For AB method we only go through once, but do so with some additional operations. 

                        if (i == 0 && iteration == 1 && validate == false && methodType == 0) {
                            //don't take unecessary steps, if we are on the first step and have no need for the large step, ignore it.
                            //Since we always want the first step to go through don't bother calculating things we don't need. 
                            iteration = 2;
                            //This doesn't actually apply to inherently adaptive methods since we cheat and do it in one iteration. 
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

                        double K[rows-methodType*quickPatch][dydt->dimension];
                        //These are the K-values that are required to evaluate RK-like methods. 
                        //They will be determined based on the provided butcher table.
                        //This is a 2D matrix since each diffyQ has its own set of K-values. 
                        //Note that we subtract the method type from the row: adaptive RK butcher tables are larger. 

                        //Since we'll be calling K while it's empty, even though there should be no errors due
                        //to the way it's set up, let's go ahead and fill it with zeroes.
                        for (int j = 0; j<rows-methodType*quickPatch; j++) {
                            for (int n = 0; n<dydt->dimension; n++) {
                                K[j][n]=0.0;
                            }
                        } 

                        double yInsert[dydt->dimension];
                        //We also need an array for the inserted y-values for each equation. 
                        //Most applications actually have the different yInsert values be independent, so 
                        //if we knew the form of the equation we could simplify the code.
                        //However, we need to make sure to always fill everything in case we have a system
                        //of the form y'=f(u,y) u'=g(u,y)

                        double dyOut[dydt->dimension];
                        //GSL demands that we use two separate arrays for y and y', so here's y'. 

                        double cInsert[numberOfConstants];

                        struct constantParameters cpInsert; 

                        updateConstantParameters(&cpInsert, &cp);
                        //Create an array to hold the constants we want.
                        //Find way to generalize. 

                        for (int j = 1; j < rows-methodType*quickPatch; j++) {
                            //Due to the way the Butcher Table is formatted, start our index at 1 and stop at the end. 
                            double xInsert = currentPosition+shift*step*scale + butcher[j-1][0]*step*scale;
                            //xInsert does not change much for different tables, just adjust the "step correction" term.
                            //xInsert is the same for every equation too.

                            for (int n = 0; n < dydt->dimension; n++) {
                                yInsert[n] = ySmolSteps[n];
                            } 
                            //Note that we are setting our buffer value, yInsert, to ySmolSteps. 
                            //This is because ySmolSteps is y at first, but we will need to evolve it forward
                            //two steps, so on the second small step this will be different. 

                            for (int n = 1; n < columns; n++) {
                                //Once again, start at index of 1 rather than 0.
                                for (int q = 0; q < dydt->dimension; q++) {
                                    yInsert[q] = yInsert[q] + butcher[j-1][n]*K[n][q];
                                    
                                }
                                //Each individual yInsert portion is dependent on one of the K values.
                                //K values are initially set to zero even though technically whenever 
                                //we would use an undeclared K-value the butcher table would have a zero.
                                //You know, just in case something goes wrong. 
                            }
                            
                            //Check for any limitations on our results. 
                            exceptionHandler(xInsert,yInsert,&cpInsert);

                            //Evaluate the constants. 
                            constEval(xInsert, yInsert,&cpInsert);

                            //Now we actually evaluate the differential equations.
                            system.function(xInsert, yInsert, dyOut, &cpInsert);
                            //yInsert goes in, dyOut comes out. Originally yInsert was overridden, this no longer happens. 

                            for (int n = 0; n < dydt->dimension; n++) {
                                K[j][n] = step*scale*dyOut[n];
                                //Fill in the K-values. 
                            } 
                        }

                        //Now that we have all the K-values set, we need to find the actual result in one final loop.
                        for (int n = 0; n< dydt->dimension; n++) {
                            K[0][n] = ySmolSteps[n]; //The 0th spot in the K-values is reserved for holding the 
                            //final value while it's being calculated. 
                            for (int j = 1; j < columns; j++) {
                                K[0][n] = K[0][n] + butcher[rows-1-methodType*quickPatch][j]*K[j][n]; 
                                //This is where the actual approximation is finally performed. 
                            }
                            ySmolSteps[n] = K[0][n]; //Set ySmol to the new estimated value. 
                        }
                        //Note that we specifically set ySmol to the value, not anything else. 
                        //This is because we wish to avoid abusing if statements and would like to do the below if only once.
                        exceptionHandler(currentPosition+(1.0+shift)*step*scale,ySmolSteps,&cpSmolSteps);
                        constEval(currentPosition+(1.0+shift)*step*scale,ySmolSteps,&cpSmolSteps); 
                        //Check for exceptions and evaluate constants. 

                        if (iteration == 1 || (i == 0 && validate == true && methodType == 1)) {
                            for (int n = 0; n<dydt->dimension; n++) {
                                yBigStep[n] = ySmolSteps[n];
                                ySmolSteps[n] = y[n];

                                //we still need to reset the value for SmolSteps on the first iteration
                                //no matter the type of method. 
                            }
                        }
                        //This only runs on the first iteration (or during validation), setting the big step to the right value
                        //and resetting the small steps for when we actually use it. 
                        //This odd structure exists purely for efficiency. 
                        
                        //If we are in an adaptive method situation, use that method and exit the iterations loop.
                        if (methodType == 1) {
                            for (int n = 0; n< dydt->dimension; n++) {
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
                                iteration = 4; //break out after we get to the end, we don't need to go any further. 
                                //We avoid this on the first iteration only if we want to validate. 
                            }
                        }

                        if (validate == true && i == 0) {
                            if (methodType == 0 && iteration == 2) {
                                //The normal validation for the normal RK methods. 

                                //Now that we've performed the approximation's first step at half the size, we can estimate the order.

                                //Create 2 arrays to hold the true values. 
                                double truthValidateBig[dydt->dimension];
                                double truthValidateSmol[dydt->dimension];
                                //Fill it with the true values. 
                                knownQEval(bound+step,truthValidateBig);
                                knownQEval(bound+step*0.5,truthValidateSmol);
                                //Then from this calculate the estimated errors.

                                for (int n = 0; n < dydt->dimension; n++) {
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
                                    for (int n = 0; n < dydt->dimension; n++) {
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
                                    double truthValidateBig[dydt->dimension];
                                    double truthValidateSmol[dydt->dimension];
                                    //Fill it with the true values. 
                                    knownQEval(bound+step,truthValidateBig);
                                    knownQEval(bound+step*0.5,truthValidateSmol);
                                    
                                    //Then from this calculate the estimated errors.
                                    for (int n = 0; n < dydt->dimension; n++) {
                                        truthValidateBig[n] = (truthValidateBig[n] - yBigStepHalf[n]);
                                        truthValidateSmol[n] = (truthValidateSmol[n] - yBigStep[n]);
                                        //Now the validation steps contain their own errors, we can compare them.
                                        printf("Order of ErrorA: %i\t%f\n",n, log2(truthValidateBig[n]/truthValidateSmol[n]));
                                        //print out the estimated error. 
                                    }
                                    knownQEval(bound+step,truthValidateBig);
                                    knownQEval(bound+step*0.5,truthValidateSmol);

                                    for (int n = 0; n < dydt->dimension; n++) {
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
                                    for (int n = 0; n < dydt->dimension; n++) {
                                        yBigStep[n] = yBigStepHalf[n];
                                        ySmolSteps[n] = ySmolStepsHalf[n];
                                        //now undo what we just did so the corret values can be read.
                                    }
                                }
                            } 

                        }

                        if (methodType == 2) {
                            iteration = 4;
                            //We only iterate once for AB. 
                            for (int n = 0; n < dydt->dimension; n++) {
                                ySmolSteps[n] = yBigStep[n];
                            }
                        }
                    }
                    //Now that the step and double step have been taken (or whatever alternative method is being used),
                    //time to calculate some errors and see if we move on to the next step. 
                    //First, from our parameters declared at the beginning, determine what our error limit is. 
                    //Using GSL's version we frist estimate our errro based on what we know.
                    if (i != 0 && methodType != 2) {
                        //Literally none of this is used for the AB method. 
                        for (int n = 0; n<dydt->dimension; n++) {
                            errorEstimate[n] = sqrt((yBigStep[n] - ySmolSteps[n])*(yBigStep[n] - ySmolSteps[n]))* errorSafety;
                            //The 4/15 is taken from GSL's solver, a 'saftey factor' with unknown reasoning. 
                        }
                        errorEstimate[1+dydt->dimension] = sqrt((cpBigStep.rho - cpSmolSteps.rho)*(cpBigStep.rho - cpSmolSteps.rho))*errorSafety;
                        //find a way to generalize to any number of constants in the struct. 

                        double errorLimiter[dydt->dimension];
                        //since the definition of the error limiter uses a derivative, we cannot use it to limit the constant's error. 
                        //We originally had the error limiter set its own values. 
                        //GSL's formatting requries us to change this. 
                        system.function(currentPosition+step,ySmolSteps, errorLimiter, &cp);
                        //Now SmolSteps is used to set the errorLimiter. 
                        for (int n = 0; n<dydt->dimension; n++) {
                            errorLimiter[n] = absoluteErrorLimit + relativeErrorLimit*(ayErrorScaler*sqrt(ySmolSteps[n]*ySmolSteps[n]) + adyErrorScaler*step*sqrt(errorLimiter[n]*errorLimiter[n]));
                        }
                        
                        //The error limiter is set for every equation. Now we need to perform checks.

                        double ratioED = 0.0;
                        for (int n = 0; n<dydt->dimension; n++) { 
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
                        if (noAdaptiveTimestep == false) {
                            //If we have no trouble...
                            if (underError == false && overError == false) {
                                errorSatisfactory = true;
                            }
                            //...say that we're cleared to move to the next step. However, if one of them was triggered, we need to adjust. 
                            //In these cases we change the actual step size. 
                            //It is theoretically possible for both to be triggered on different equations. In that case, overError
                            //takes prescedent. We would rather have more accuracy than less in odd situations like that. 

                            //These if statements perform step adjustment if needed. Based on GSL's algorithm. 
                        
                            else if (overError == true) {
                                step = step * scaleFactor * pow(ratioED,-1.0/butcher[rows-1-methodType*quickPatch][0]);
                            } else { //if underError is true and overError is false is the only way to get here. The true-true situation is skipped.
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
                        //With that, the step size has been changed. If errorSatisfactory should be false, it goes back and performs everything again
                        //with the new step size. 
                    } else {
                        errorSatisfactory = true;
                        //We always want the *first* step to go through without change, often the first step is chosen for a specific reason. 
                        //In our work this generally came from a need to plot data sets against each other. 
                        //Also do this if we are using the AB method, as it has no error checks. 
                    }
                }
                
                //Finally, we actually update the real answer. 
                for (int n = 0; n<dydt->dimension; n++) {
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
                for (int n = 0; n < dydt->dimension; n++) {
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
                    for (int n = 0; n < dydt->dimension; n++) {
                        fprintf(fp, "Equation %i:,\t%15.14e,\t",n,errorEstimate[n]);
                    }
                    for (int n = 0; n < numberOfConstants; n++) {
                        fprintf(fp, "Constant %i:,\t%15.14e,\t",n,errorEstimate[n+dydt->dimension]);
                    }   
                    fprintf(fp,"\n");
                }
                
                if (reportErrorActual == true) {
                    //Now if we have an actual error to compare against with, there's some more work to do. 
                    double yTruth[dydt->dimension];
                    double cTruth[numberOfConstants];
                    struct constantParameters cpTruth; 
                    //True values for everything we compare with.

                    updateConstantParameters(&cpTruth, &cp);
                    //Setting of constants (generalize?)

                    knownQEval(currentPosition,yTruth);
                    constEval(currentPosition,yTruth,&cpTruth);
                    assignConstants(c,&cp); 
                    assignConstants(cTruth,&cpTruth); 
                    fprintf(fp, "Errors:,\t");
                    for (int n = 0; n < dydt->dimension; n++) {
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

                if (stepType == nrpy_odiegm_step_AB) {
                    //At the END of every loop, we "shift" the values in the array "down" one space, that is, into the "past"
                    //Present values are 0, previous step is 1, step before that is 2, etc. 
                    for (int n = 0; n < dydt->dimension; n++) {
                        for (int m = adamsBashforthOrder - 1; m > 0; m--) {
                            yValues[n][m] = yValues[n][m-1];
                            //note that we start at the last column, m, and move the adjacent column to it. 
                            //This pushes off the value at the largest m value, since it's far enough in the past we no longer care.
                        }
                        yValues[n][0] = y[n]; 
                        //present values update to what we just calculated. 
                        //We have now completed stepping. 
                    }  
                }

                //And the very last thing we do in the loop is ask if we terminate it. 
                if (doWeTerminate(currentPosition, y, &cp) == 1) {
                    i = SIZE;
                }
    //This will presumably take forever to actually complete. 
    return 0;                      
}
int nrpy_odiegm_evolve_apply_fixed_step (nrpy_odiegm_evolve * e,
                                        nrpy_odiegm_control * con,
                                        nrpy_odiegm_step * step,
                                        const nrpy_odiegm_system * dydt,
                                        double *t, const double h0,
                                        double y[]){
    //This method performs a ton of time steps at the same step size. 

    return 0;
}
int nrpy_odiegm_driver_apply (nrpy_odiegm_driver * d, double *t,
                             const double t1, double y[]){
    //Takes a step at the driver level. 
    //Seems unecessary to us.

    return 0;
}
int nrpy_odiegm_driver_apply_fixed_step (nrpy_odiegm_driver * d, double *t,
                                        const double h,
                                        const unsigned long int n,
                                        double y[]){
    //Takes several steps of identical size at the driver level.  
    //Seems unecessary to us.

    return 0;
}

