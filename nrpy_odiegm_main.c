#include "nrpy_odiegm_funcs.c" // nrpy_odiegm itself.
#include "nrpy_odiegm_user_methods.c" // User-dependent functions. 

// This file is technically not part of Odie, it is just an example implementation.
// However, it is exceptionally versatile, and can be used to run virtually 
// every system of differential equations desired.
// That said, since it's general, it trades off some efficiency to handle
// every type of method. If high efficiency is desired, it is recommended that the user 
// make a custom main function. 

int main()
{
    printf("Beginning ODE Solver \"Odie\" V10...\n");

    // SECTION I: Preliminaries

    // Before the program actually starts, variables need to be created
    // and set, as well as the functions chosen. 
    // The system of differential equations can be found declared in diffy_Q_eval
    // in nrpy_odiegm_user_methods.c

    double step = 0.00001; // the "step" value. Initial step if using an adaptive method.
    double current_position = 0.0; // where the boundary/initial condition is. 
    // Same for every equation in the system.
    int number_of_equations = 4; // How many equations are in our system?
    int number_of_constants = 1; // How many constants do we wish to separately evaluate and report? 
    // If altering the two "numberOf" ints, be careful it doesn't go over the actual number 
    // and cause an overflow in the functions in nrpy_odiegm_user_methods.c
    const int size = 100000; // How many steps are we going to take? 
    // This is the default termination condition. 
    int adams_bashforth_order = 4; // If using the AB method, specify which order you want.
    // If we are not using the AB method this is set to 0 later automatically. 4 by default. 
    bool no_adaptive_step = false; // Sometimes we just want to step forward uniformly 
    // without using GSL's awkward setup. False by default. 

    bool report_error_actual = false;
    bool report_error_estimates = false;
    // AB methods do not report error estimates. 
    // BE WARNED: setting reporError (either kind) to true makes
    // it print out all error data on another line,
    // the file will have to be read differently. 

    // ERROR PARAMETERS: Use these to set limits on the erorr. 
    double absolute_error_limit = 1e-14; // How big do we let the absolute error be?
    double relative_error_limit = 1e-14; // How big do we let the relative error be?
    // Default: 1e-14 for both.
    // Note: there are a lot more error control numbers that can be set inside the 
    // control "object" (struct) d->c.

    char file_name[] = "ooData.txt"; // Where do you want the data to print?

    // Now we set up the method. 
    const nrpy_odiegm_step_type * step_type;
    step_type = nrpy_odiegm_step_AB;
    // Here is where the method is actually set, by specific name since that's what GSL does. 

    const nrpy_odiegm_step_type * step_type_2;
    step_type_2 = nrpy_odiegm_step_AB;
    // This is a second step type "object" (struct) for hybridizing. 
    // Only used if the original type is AB.
    // Set to AB to use pure AB method. 

    // AFTER THIS POINT THERE SHOULD BE NO NEED FOR USER INPUT, THE CODE SHOULD HANDLE ITSELF. 

    // We need to define a struct that can hold all possible constants. 
    struct constant_parameters cp; 
    cp.dimension = number_of_constants;
    // We'll set the actual parameters later. 
    // Do note that cp itself needs to be declared in constant_parameters in 
    // nrpy_odiegm_user_methods.c manually.
    // The methods that make use of it it need to be declared as well, if they are used.

    nrpy_odiegm_system system = {diffy_Q_eval,known_Q_eval,number_of_equations,&cp};
    // This is the system of equations we solve.
    // The second slot was originally the Jacobian in GSL, but we use it to pass a 
    // true answer function that may or may not be used.

    nrpy_odiegm_driver *d;
    d = nrpy_odiegm_driver_alloc_y_new(&system, step_type, step, absolute_error_limit, relative_error_limit); 
    // This is the "object" (struct) that runs everything, contains every needed varaible, etc. 
    // Basically the master of the whole thing, hence why it's called the "driver"
    // Contains three major sub-objects besides the step type. 
    // c is the controller, which is primarily used to store adaptive timestep values. 
    // s is the step, which has the step type in it, but also parameters that describe the steps.
    // e is the evolver, which actually performs the update when it is requested. 

    int method_type = 1;
    if (step_type->rows == step_type->columns) {
        method_type = 0; // AKA, normal RK-type method. 
    } // No need for an else, we set it to 1 earlier to represent Adaptive methods. 
    if (step_type->rows == 19) { 
        method_type = 2;
    } else {
        adams_bashforth_order = 0;
    }
    d->s->adams_bashforth_order = adams_bashforth_order;
    d->e->no_adaptive_step = no_adaptive_step;
    // Based on what type of method we are using, we adjust some parameters within the driver.

    if (method_type == 2) {
        printf("Method Order: %i.\n",adams_bashforth_order);
    } else {
        printf("Method Order: %i.\n",step_type->order);            
    }
    
    double y[number_of_equations];
    // These next few variables temporarily store the values calculated before they are 
    // printed to the output file and forgotten.
    // y contains the values of the actual equations. 
    // Each array only holds values at one evaluation point, but one for each Equation.

    double c[number_of_constants];
    // c is just used to hold any constants we wish to report. 
    // You'd think that, since we have the constants in a struct, we can avoid declaring this.
    // No. Not as far as we can tell, anyway. Structs are a pain to iterate through,
    // and we can't know what form the user is going to hand us the struct in. 

    // This here sets the initial conditions as declared in get_initial_condition
    get_initial_condition(y); 
    const_eval(current_position, y,&cp);

    FILE *fp2;
    fp2 = fopen(file_name,"w");
    printf("Printing to file '%s'.\n",file_name);

    // Open the file we'll be writing data to. 
    // Before continuing, let's print out our initial data. 
    // The print function is automatically adaptable to any size of data. 
    // We print both to the terminal and to the file for the initial conditions, 
    // but later only print to the file.

    // First, print the location we are at. 
    printf("INITIAL: Position:,\t%f,\t",current_position);
    fprintf(fp2, "Position:,\t%15.14e,\t",current_position);
    // Second, go through and print the result for every single equation in our system.
    for (int n = 0; n < number_of_equations; n++) {
        printf("Equation %i:,\t%15.14e,\t",n, y[n]);
        fprintf(fp2, "Equation %i:,\t%15.14e,\t",n, y[n]);
    }
    // Third, print out desired constants.
    assign_constants(c,&cp);    
    for (int n = 0; n < number_of_constants; n++) {
        printf("Constant %i:,\t%15.14e,\t",n, c[n]);
        fprintf(fp2, "Constant %i:,\t%15.14e,\t",n, c[n]);
    }
    // Lastly, the newline character. 
    printf("\n");
    fprintf(fp2,"\n");
    // Comma delimiters are printed to the file so it can be read as .csv with ease. 

    if (report_error_estimates == true) {
        // In order to keep things neat and regular in the file, print a first line of errors. 
        // Even though by necessity all of them must be zero. 
        fprintf(fp2, "Errors Estimates:,\t");
        for (int n = 0; n < number_of_equations; n++) {
            fprintf(fp2, "Equation %i:,\t0.0,\t",n);
        }
        for (int n = 0; n < number_of_constants; n++) {
            fprintf(fp2, "Constant %i:,\t0.0,\t",n);
        }   
        fprintf(fp2,"\n");
    }
    
    if (report_error_actual == true) {
        // In order to keep things neat and regular in the file, print a first line of errors. 
        // Even though by necessity all of them must be zero. 
        fprintf(fp2, "Errors:,\t");
        for (int n = 0; n < number_of_equations; n++) {
            fprintf(fp2, "Equation %i:,\t0.0,\t",n);
            fprintf(fp2, "Truth:,\t%15.14e,\t",y[n]);
        }
        for (int n = 0; n < number_of_constants; n++) {
            fprintf(fp2, "Constant %i:,\t0.0,\t",n);
            fprintf(fp2, "Truth:,\t%15.14e,\t",c[n]);
        }   
        fprintf(fp2,"\n");
    }

    // SECTION II: The Loop

    // This loop fills out all the data.
    // It takes a provided butcher table and executes the method stored within. 
    // Any RK table should work, even one not included by default.
    // Also handles AB methods up to 19th order. No one should ever need more. 

    for (int i = 0; i < size; i++){
        
        // Hybrid Methods require some fancy footwork, hence the if statements below. 
        if (method_type == 2 && i == 0 && step_type_2 != nrpy_odiegm_step_AB) {
            d->s->type = step_type_2;
            d->s->rows = step_type_2->rows;
            d->s->columns = step_type_2->columns;
            d->s->method_type = 0;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        } else if (step_type != step_type_2 && method_type == 2 && i == adams_bashforth_order) {
            d->s->type = step_type;
            d->s->rows = step_type->rows;
            d->s->columns = step_type->columns;
            d->s->method_type = 2;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        }

        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, &system, &current_position, current_position+step, &step, y);
        // This is the line that actually performs the step.

        exception_handler(current_position,y);
        const_eval(current_position,y,&cp);
        assign_constants(c,&cp);
        // These lines are to make sure the constant updates. 
        // And exception constraints are applied.  

        // Printing section.
        // Uncomment for live updates. Prints to the file automatically.
        // printf("Position:,\t%15.14e,\t",current_position);
        fprintf(fp2, "Position:,\t%15.14e,\t",current_position);
        for (int n = 0; n < number_of_equations; n++) {
            // printf("Equation %i:,\t%15.14e,\t",n, y[n]);
            fprintf(fp2, "Equation %i:,\t%15.14e,\t",n, y[n]);
        }

        for (int n = 0; n < number_of_constants; n++) {
            // printf("Constant %i:,\t%15.14e,\t",n, c[n]);
            fprintf(fp2, "Constant %i:,\t%15.14e,\t",n, c[n]);
            // printf("Constant %i:,\t%15.14e %15.14e,\n",n, c[n], y[n]);
        }
        // printf("\n");
        fprintf(fp2,"\n");

        if (report_error_estimates == true) {
            // Print the error estimates we already have. 
            fprintf(fp2, "Error Estimates:,\t");
            for (int n = 0; n < number_of_equations; n++) {
                fprintf(fp2, "Equation %i:,\t%15.14e,\t",n,*(d->e->yerr)); 
            }
            // Constant estimates not reported, only differential equation values. 
            fprintf(fp2,"\n");
        }
            
        if (report_error_actual == true) {
            // Now if we have an actual error to compare against, there's some more work to do. 
            double y_truth[number_of_equations];
            double c_truth[number_of_constants];
            struct constant_parameters cp_truth; 
            // True values for everything we compare with.
            
            known_Q_eval(current_position,y_truth);
            const_eval(current_position,y_truth,&cp_truth);

            assign_constants(c,&cp); 
            assign_constants(c_truth,&cp_truth);
 
            fprintf(fp2, "Errors:,\t");
            for (int n = 0; n < number_of_equations; n++) {
                fprintf(fp2, "Equation %i:,\t%15.14e,\t",n, y_truth[n]-y[n]);
                fprintf(fp2, "Truth:,\t%15.14e,\t",y_truth[n]);
            }
            for (int n = 0; n < number_of_constants; n++) {
                fprintf(fp2, "Constant %i Error:,\t%15.14e,\t",n, c_truth[n]-c[n]);
                fprintf(fp2, "Truth:,\t%15.14e,\t",c_truth[n]);
            } 
            fprintf(fp2,"\n");
        }

        if (do_we_terminate(current_position, y, &cp) == 1) {
            i = size;
            // If we need to bail, set i to size to break the loop.
            // Also potentially a good idea: print the final line. 
            printf("FINAL: Position:,\t%15.14e,\t",current_position);
            for (int n = 0; n < number_of_equations; n++) {
                // printf("Equation %i:,\t%15.14e,\t",n, y[n]);
                printf("Equation %i:,\t%15.14e,\t",n, y[n]);
            }

            for (int n = 0; n < number_of_constants; n++) {
                // printf("Constant %i:,\t%15.14e,\t",n, c[n]);
                printf("Constant %i:,\t%15.14e,\t",n, c[n]);
                // printf("Constant %i:,\t%15.14e %15.14e,\n",n, c[n], y[n]);
            }
            // printf("\n");
            printf("\n");
        }
    }

    // SECTION III: Analysis

    // Minor post-processing goes here. 
    // Anything advanced will need to be done in a data analysis program. 
    // We like to use matplotlib for python.

    fclose(fp2);

    nrpy_odiegm_driver_free(d);
    // MEMORY SHENANIGANS

    printf("ODE Solver \"Odie\" V10 Shutting Down...\n");
    return 0;

}

// - GM, master of dogs.
