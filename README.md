# ODESolver
This repository houses an Ordinary Differential Equations solver known as "Odie" as well as a few of its applications. 

It is best for anyone who wants to know how to use it to rely on the tutorial jupyter notebook,
`NRPy+_OdieGM_Tutorial.ipynb`. All the information you could possibly want to know is found there. 

If you just want to run the C-code on its own, `nrpy_odiegm_main.c` is the file you want to run. 

By default it solves the TOV equations with an adaptive fourth order Runge-Kutta algorithm. 

The `OldFilesWithValidation` folder just holds some depricated code that I might want to refer to since it has a feature the current code no longer needs (self validation of method order), but still could be useful to know in the future.

The `TOVOdieGM` folder contains a Thorn designed for the Einstein Toolkit based on Odie that can solve the TOV Equations. This is intended for use in setting up initial data for neutron stars. This also serves as an example of a more advanced and involved application of Odie. 

-GM, master of cats.
