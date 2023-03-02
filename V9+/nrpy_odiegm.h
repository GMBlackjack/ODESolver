#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

double zeroArray[19][19] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
//The existence of this array is evidence of us throwing up our hands and plugging the tiny hole with an entire island. 
//No more memory leaks or overwrites, though! It solves the problem! 

typedef struct {
    int (*function) (double x, double y[], double dydx[], void *params);
    //The function passed to this struct contains the definitions of the differnetial equations. 
    //int (*jacobian) (double t, const double y[], double *dfdy, double dfdt[], void *params); //The jacobian is a holdover from GSL, it will not be used in this program,
                   //Always pass NULL to this part. 
                   //If we were not making a drop-in replacement, we would just have a three-paraemeter struct.
    int (*trueFunction) (double x, double y[]);
    //INSTEAD we will use that slot to allow a passing of a true value! 
    size_t dimension; //For storing how big our system of equations is. 
    //Just pass it an int, usually. 
    void *params; //For storing extra constants needed to evaluate the functions. 
    //params->dimension stores how many there are. 
} nrpy_odiegm_system;

typedef struct {
    //unlike with the system object above, this step_type object does not need
    //to match GSL's form explicitly, it just needs to define the method.
    int rows; 
    int columns; //since we are passing a void pointer to do this, we need a way
    int order; //record the order, stop with the sillyness if we can. 
    //to know how large it is in the end.
    //These are set in butcher.c.
    //void *(*alloc) (size_t dim); //Adding this breaks the array, we can't have it. 
    void *butcher;
    //Make sure to put this at the end of the struct though,
    //in case we add more parts to it, nonspecific arrays must be the last element.

    //two of these step_type objects might be needed at once, depending on implementation. 
} nrpy_odiegm_step_type;

typedef struct {
  const nrpy_odiegm_step_type *type; 
  int rows; 
  int columns; //since we are passing a void pointer to do this, we need a way
  //to know how large it is in the end.
  int methodType; //can be calculated from above, but is easier to just store. 0,1,2 values. 
  int adamsBashforthOrder; //order if an AB method is used.
  void *yValues; //The extremely funky parameter that hides a 2D array, used when
  //the past steps are important for AB method.  
  //Stored in step since it needs access to adamsBashforthOrder for allocation.
} nrpy_odiegm_step;

typedef struct {
    double absLim;
    double relLim;
    double scaleFactor;
    double errorSafety;
    double ayErrorScaler;
    double adyErrorScaler;
    double maxStepAdjustment;
    double minStepAdjustment;
    double absoluteMaxStep;
    double absoluteMinStep;
    double errorUpperTolerance;
    double errorLowerTolerance;
    //We added these ourselves. Control the error!
    //we suppose this means that our control object acts NOTHING like GSL's control object
    //save that it stores error limits. 
} nrpy_odiegm_control;

typedef struct
{
  double *y0; 
  double *yerr; 
  //double *dydt_in;
  //double *dydt_out; //These two values are not used as we are not using Jacobian methods. 
  double last_step;
  double bound; //the point at which we started is sometimes important. 
  double currentPosition;
  unsigned long int count; //equivalent to i. 
  //unsigned long int failed_steps; //we automatically handle tossing out "bad" steps.
  bool noAdaptiveTimestep;
  bool validate;
} nrpy_odiegm_evolve;

typedef struct {
    const nrpy_odiegm_system *sys; /* ODE system */
    nrpy_odiegm_evolve *e;         /* evolve object */
    nrpy_odiegm_control *c;
    nrpy_odiegm_step *s;
    double h;                     /* step size */
    //Curiously, this is where the step size is held. Usually it's passed to functions directly though. 
} nrpy_odiegm_driver;

//Allocation Functions below. 

nrpy_odiegm_step *
nrpy_odiegm_step_alloc (const nrpy_odiegm_step_type * T, size_t dim)
{
  nrpy_odiegm_step *s = (nrpy_odiegm_step *) malloc (sizeof (nrpy_odiegm_step));
  s->type = T;
  s->methodType = 1;
  s->adamsBashforthOrder = 0;
  if (T->rows == T->columns) {
    s->methodType = 0; //aka, normal RK-type method. 
  }
  if (T->rows == 19) {
    s->methodType = 2; //AB method. 
    s->adamsBashforthOrder = 4; //default order chosen, if user wants control they will specify elsewhere. 
  }

  double emptyArray[dim][s->adamsBashforthOrder];
  for (int n = 0; n < dim; n++) {
    for (int m = 0; m < s->adamsBashforthOrder; m++) {
      emptyArray[n][m] = 0;
    }
  }
  s->yValues = &zeroArray;

  return s;
}

nrpy_odiegm_evolve *
nrpy_odiegm_evolve_alloc (size_t dim)
{
  nrpy_odiegm_evolve *e = (nrpy_odiegm_evolve *) malloc (sizeof (nrpy_odiegm_evolve));
  e->y0 = (double *) malloc (dim * sizeof (double));
  e->yerr = (double *) malloc (dim * sizeof (double));
  //e->dimension = dim; //this value turns out to be completely unecessary as the system stores the dimension. 
  e->count = 0;

  e->last_step = 0.0; 
  e->bound = 0.0; //This will need to be adjusted to handle other starting positions. 
  e->currentPosition = e->bound;
  e->noAdaptiveTimestep = false; //for now, this one will be set with the other methods. 
  e->validate = false; //this one will need to be manually set 
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
    //These are all the default values.

  return c;
}

nrpy_odiegm_driver * nrpy_odiegm_driver_alloc_y_new (const nrpy_odiegm_system * sys,
                               const nrpy_odiegm_step_type * T,
                               const double hstart,
                               const double epsabs, const double epsrel)
{
    //Initializes an ODE driver system with control object of type y_new.

    nrpy_odiegm_driver *state;
    state = (nrpy_odiegm_driver *) calloc (1, sizeof (nrpy_odiegm_driver)); //valgrind doesn't like this line. 
    const size_t dim = sys->dimension; 
    state->sys = sys;
    state->s = nrpy_odiegm_step_alloc (T, dim);

    state->e = nrpy_odiegm_evolve_alloc (dim);
    state->h = hstart;

    state->c = nrpy_odiegm_control_y_new (epsabs, epsrel);
    //Otherwise we're not using a method that needs step size control. 
  //There were functions here that assigned the driver to the objects contained in the driver.
  //We will not be doing that nonsense. 

  return state;
}

//Memory freeing methods. 
//Only necessary since we want a drop in replacement.
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