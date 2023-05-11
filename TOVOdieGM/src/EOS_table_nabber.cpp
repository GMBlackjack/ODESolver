// This file has one purpose and one purpose only: grab EOS_Omni information and tabulate it
// So the rest of Odie can read it. 
// the first part here is essentially copied from the EOS_Omni header
// verbatum, just with unecessary parts removed.

#include <iostream>
#include <cmath>
#define NUC_EOS_HH

#include "cctk.h"

#define HAVEGR 1
#define NTABLES 19

#define mu_nu_index 3
#define logpress_index 0
#define logenergy_index 1

namespace nuc_eos_private {

// table data

  extern int nrho;
  extern int ntemp;
  extern int nye;

  extern double * restrict alltables;
  // extern double * restrict epstable;
  extern double * restrict logrho;
  extern double * restrict logtemp;
  /* extern double dlintemp,dlintempi;
  extern double drholintempi;
  extern double dlintempyei;
  extern double drholintempyei;
  extern double * restrict yes;
  extern double dtemp, dtempi;
  extern double drho, drhoi;
  extern double dye, dyei;
  extern double drhotempi;
  extern double drhoyei;
  extern double dtempyei;
  extern double drhotempyei; */

}

using namespace std;
using namespace nuc_eos_private;

extern "C" void EOS_table_values_nabber(double T_initial, double **newlogrho, double **newlogpres, double **newlogeps, int *array_size)
{
    cout.precision(16); // Make sure we have the right number of sigs. 
    
    // first find out the closest index for T. 
    double multiUseBuffer = fabs(log(T_initial) - logtemp[0]);
    int tempIndex = 0;
    for (int i = 1; i < ntemp-1;i++) {
    	if (multiUseBuffer > fabs(log(T_initial) - logtemp[i])) {
    		tempIndex = i;
    		multiUseBuffer = fabs(log(T_initial) - logtemp[i]);
    	}
    }
    
    // Also find what we need for the electron fraction index.
    // Admittedly did this by following suggestions by Leo Werneck.
    // Don't really undersatnd why it's all arranged like this. 
    int electronFracIndex[nrho];
    for (int i = 0; i < nrho-1;i++) {
    	multiUseBuffer = alltables[mu_nu_index + NTABLES*(i + nrho*((tempIndex) + ntemp*(0)))];
    	electronFracIndex[i] = 0;
    	for(int k = 1; k < nye-1;k++) {
    		if (fabs(multiUseBuffer) > fabs(alltables[mu_nu_index + NTABLES*(i + nrho*((tempIndex) + ntemp*(k)))])) {
    			electronFracIndex[i] = k;
    			multiUseBuffer = alltables[mu_nu_index + NTABLES*(i + nrho*((tempIndex) + ntemp*(k)))];
    		}
    	}
    }
    
    // Now we've chosen the part of the table we want to look at. 
    // And so we now set up the actual arrays we want to read. 

    double *rho_local = (double*)malloc(nrho * sizeof(double));
    double *pres_local = (double*)malloc(nrho * sizeof(double));
    double *eps_local = (double*)malloc(nrho * sizeof(double));
    for (int i = 0; i < nrho-1;i++) {
    	rho_local[i] = logrho[i];
    	pres_local[i] = alltables[logpress_index + NTABLES*(i + nrho*((tempIndex) + ntemp*(electronFracIndex[i])))];
    	eps_local[i] = alltables[logenergy_index + NTABLES*(i + nrho*((tempIndex) + ntemp*(electronFracIndex[i])))];
    }
    
    // Set the arrays to these ones so the C-code can access them.
    // Remember, they need to be freed by the C-code later. 
    *newlogrho = rho_local;
    *newlogpres = pres_local;
    *newlogeps = eps_local;
    *array_size = nrho; // We also need to know how big to make the arrays. 
    
}
