#include "nrpy_odiegm_funcs.c" //nrpy_odiegm itself.
#include "nrpy_odiegm_user_methods.c" //user-dependent functions. 

int main()
{

    double step = 0.01; 
    double bound = 0.0;
    int numberOfEquations = 4; 
    int numberOfConstants = 1; 
    const int SIZE = 10000; 
    double absoluteErrorLimit = 1e-14; 
    double relativeErrorLimit = 1e-14; 

    const nrpy_odiegm_step_type * stepType;
    stepType = nrpy_odiegm_step_RK4;
    struct constantParameters cp; 
    cp.dimension = numberOfConstants;
    nrpy_odiegm_system system = {diffyQEval,knownQEval,numberOfEquations,&cp};

    nrpy_odiegm_driver *d;
    d = nrpy_odiegm_driver_alloc_y_new(&system, stepType, step, absoluteErrorLimit, relativeErrorLimit); 
    d->e->bound = bound;

    double y[numberOfEquations];
    double c[numberOfConstants];
    double currentPosition = bound;

    getInitialCondition(y); 
    constEval(bound, y,&cp);

    for (int i = 0; i < SIZE; i++){
        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, &system, &currentPosition, currentPosition+step, &step, y);
        exceptionHandler(currentPosition,y);
        constEval(currentPosition,y,&cp);
        assignConstants(c,&cp);

        printf("Position:,\t%15.14e,\t",currentPosition);
        for (int n = 0; n < numberOfEquations; n++) {
            printf("Equation %i:,\t%15.14e,\t",n, y[n]);
        }
        for (int n = 0; n < numberOfConstants; n++) {
            printf("Constant %i:,\t%15.14e,\t",n, c[n]);
        }
        printf("\n");
        if (doWeTerminate(currentPosition, y, &cp) == 1) {
            i = SIZE;
        }
    }

    nrpy_odiegm_driver_free(d);
    return 0;
}

// - GM, master of dogs.