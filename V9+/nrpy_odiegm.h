#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

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
} nrpy_odiegm_system;

typedef struct {
    //unlike with the system object above, this step_type object does not need
    //to match GSL's form explicitly, it just needs to define the method.
    //Since we are using RK-type methods exclusively, this is just a
    int rows; 
    int columns; //since we are passing a void pointer to do this, we need a way
    //to know how large it is in the end.
    void *butcher;
    //Make sure to put this at the end of the struct though,
    //in case we add more parts to it, nonspecific arrays must be the last element.
} nrpy_odiegm_step_type;

typedef struct {
  //const gsl_odeiv2_step_type *type;
  size_t dimension;
  void *state;
} nrpy_odiegm_step;

typedef struct {
  //const gsl_odeiv2_control_type *type;
  void *state;
} nrpy_odiegm_control;

typedef struct
{
  size_t dimension;
  double *y0;
  double *yerr;
  double *dydt_in;
  double *dydt_out;
  double last_step;
  unsigned long int count;
  unsigned long int failed_steps;
  //const nrpy_odiegm_driver *driver;
} nrpy_odiegm_evolve;

typedef struct {
    const nrpy_odiegm_system *sys; /* ODE system */
    nrpy_odiegm_step *s;           /* stepper object */
    nrpy_odiegm_control *c;        /* control object */
    nrpy_odiegm_evolve *e;         /* evolve object */
    double h;                     /* step size */
    double hmin;                  /* minimum step size allowed */
    double hmax;                  /* maximum step size allowed */
    unsigned long int n;          /* number of steps taken */
    unsigned long int nmax;       /* Maximum number of steps allowed */
} nrpy_odiegm_driver;

//High-level functions. 

nrpy_odiegm_driver * nrpy_odiegm_driver_alloc_y_new (const nrpy_odiegm_system * sys,
                               const nrpy_odiegm_step_type * T,
                               const double hstart,
                               const double epsabs, const double epsrel)
{
    /* Initializes an ODE driver system with control object of type y_new. */

    nrpy_odiegm_driver *state;

    state = (nrpy_odiegm_driver *) calloc (1, sizeof (nrpy_odiegm_driver));

    const size_t dim = sys->dimension;

    if (dim == 0)
    {
        //gsl_odeiv2_driver_free(state);
    }

    state->sys = sys;

    state->s = nrpy_odiegm_step_alloc (T, dim);

    if (state->s == NULL)
      {
        //gsl_odeiv2_driver_free(state);
      }

    state->e = gsl_odeiv2_evolve_alloc (dim);
  
  if (state->e == NULL)
    {
      //gsl_odeiv2_driver_free(state);
    }

    state->h = hstart;
    state->hmin = 0.0;
    state->hmax = 1.0;
    state->nmax = 0;
    state->n = 0;
    state->c = NULL;

  if (epsabs >= 0.0 && epsrel >= 0.0)
    {
      state->c = nrpy_odiegm_control_y_new (epsabs, epsrel);

      /* if (state->c == NULL)
        {
          nrpy_odiegm_driver_free (state);
          GSL_ERROR_NULL ("failed to allocate control object", GSL_ENOMEM);
        }*/
    }
  else
    {
      //gsl_odeiv2_driver_free (state);
    }

  /* Distribute pointer to driver object */


  //these things here can probably be done manually, no need for "nonsense."
  nrpy_odiegm_step_set_driver (state->s, state);
  nrpy_odiegm_evolve_set_driver (state->e, state);
  nrpy_odiegm_control_set_driver (state->c, state);

  return state;
}

