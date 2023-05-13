# ODESolver
This repository houses an Ordinary Differential Equations solver known as "Odie" as well as a few of its applications. 

It is best for anyone who wants to know how to use it to rely on the tutorial jupyter notebooks. I recommend starting with `NRPy+_OdieGM_Quickstart.ipynb` and then going from there, as it links to the other notebooks. `NRPy+_OdieGM_Code_Regeneration.ipynb` can be used to restore the C code in its original form if it is lost for some reason. 

If you just want to run the C-code on its own, `nrpy_odiegm_main.c` is the file you want to run. By default it solves the TOV equations with an adaptive fourth order Runge-Kutta algorithm. 

The `OldFilesWithValidation` folder holds some depricated code that I might want to refer to since it has a feature the current code no longer needs (self validation of method order), but still could be useful to know in the future. It also contains some old jupyter notebooks, one of which is the "master tutorial" which had all the information in one place and was absolutely tremendous; it is kept here in case something was lost in splitting it up in four. The other notebooks are testing playgrounds that might have code I wish to refer to later. 

The `TOVOdieGM` folder contains a Thorn designed for the Einstein Toolkit based on Odie that can solve the TOV Equations. This is intended for use in setting up initial data for neutron stars. This also serves as an example of a more advanced and involved application of Odie, though not one intended to serve as a template for other uses. 

There are also a bunch of .py files in here that the notebooks rely on. These are from NRPy+ and help the notebooks run the C code. 

-GM, master of cats.
