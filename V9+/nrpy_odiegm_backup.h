#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

typedef struct{ 
    int dimension; //number that says how many we have. 
    double rho;
    //add more as necessary.  
} constantParameters;

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

