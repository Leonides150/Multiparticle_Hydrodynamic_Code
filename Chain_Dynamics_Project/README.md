# Version to study particle dynamics under symmetric and antisymmetric hydrodynamic components

## Summary
The original code was modified for studying the particle-chain dynamics in the presence of both symmetric forces e.g., via Hele--Shaw (H.S.) Quadrupoles, and swapped trajectory mechanisms and antisymmetric forces e.g., Hele--Shaw Dipoles, and near-field exluded volume potentials.

## Code Structure

The main body of the code resides in the ```main.f90``` Fortran file. The ```Makefile``` is used for compilation and need to be modified to reflect proper paths to the libraries being accessed. Finally, the ```runprog*.bash``` file is used to run an instance of the code. There are many variants of this bash script aimed at different types of runs.

##List of code variants:
1. Consists of modules to simulate interactions between deformable drops in Poiseuille flow with variable deformability (equivalent to Capillary numbers).
2. Consists of augmentation made to the framework for including wall reflection to simulate an experimental setup in collaboration with [Siva Vanapalli Lab](https://www.linkedin.com/in/siva-vanapalli-7071a848/) 
