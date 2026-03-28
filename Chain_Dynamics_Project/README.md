# Version to study particle dynamics under symmetric and antisymmetric hydrodynamic components

## Summary
The original code was modified for studying the particle-chain dynamics in the presence of both symmetric forces e.g., via Hele--Shaw (H.S.) Quadrupoles, and swapped trajectory mechanisms and antisymmetric forces e.g., Hele--Shaw Dipoles, and near-field exluded volume potentials.

## Code Structure

The main body of the code resides in the ```main.f90``` Fortran file. The ```Makefile``` is used for compilation and need to be modified to reflect proper paths to the libraries being accessed. Finally, the ```runprog*.bash``` file is used to run an instance of the code. There are many variants of this bash script aimed at different types of runs.

##List of code variants:
1. Consists of modules to simulate interactions between deformable drops in Poiseuille flow with variable deformability (equivalent to Capillary numbers). The framework could be used to simulate deformable particles in Poiseuille flow or particles interacting via symmetric (Quadrupoles say) and anti-symmetric (Dipoles say) forces. 
2. Consists of augmentation made to the framework to account for **particle-side wall interactions via the reflection criterion** to simulate an experimental setup in collaboration with [Siva Vanapalli Lab](https://www.linkedin.com/in/siva-vanapalli-7071a848/). The general setup included a thin width channel with a polymer suspension popagating under Poiseuille flow conditions.
3. Consists of code augmented with **Gaussian velocity fluctuation** to account for polydispersity of droplets in the system described above in pt. 2, a **forward-pointing short-range potential** in the front to mimic propagating front of droplets at the leading end, and a **long-range potential at the trailing end** to replicate the lagging end of the droplet suspension.
