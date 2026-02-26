# IsochroneSymplecticIntegrator_QuickExample
The program provides a simple demonstration of how to use the isochrone propagator within a symplectic integrator, offering an intuitive and pedagogical introduction to the [isochrone splitting](https://ui.adsabs.harvard.edu/abs/2025A%26A...697A.106B/abstract). As an illustrative example, it computes the motion of a test particle evolving in a Plummer potential. The symplectic scheme implemented is the ‘Leapfrog’ (a second-order symplectic scheme commonly used in Hamiltonian dynamics). 

The code can be readily adapted for projects requiring integration into more complex potentials or the use of higher-order symplectic integrators.

## Key File
- `IsochronePropagator.f90` - Implements the analytical isochrone propagator.

## Installation 
Make sure your computer has the `gfortran` compiler installed.
Then, download this project and run `make` in the directory containing the `Makefile` to compile it.

No additional packages are needed.

# Basic usage
From the `in` directory, `param.txt` lets you configure the potential parameters, integration settings, isochrone splitting parameters, and initial conditions.

Once the program is compiled, run it with: 

`./run_IsochroneSplitting_QuickExample.x in/param.txt out/YourOutputFile.dat` 

from the directory where the executable is located. The program will save the results in `out/YourOutputFile.dat`.

# Contributing
Feel free to send me an email if you have any ideas to improve the speed or reduce the numerical errors of the isochrone analytical propagator.

If you encounter a bug, please let me know.

# License
Distributed under the **MIT License**. See `LICENSE.txt` for more details.

# Attribution
This example is meant for educational purposes, but you are free to reuse, modify, or incorporate it into your own projects, provided that proper credit is given. Please cite the following DOI: [10.1051/0004-6361/202553886](https://ui.adsabs.harvard.edu/abs/2025A%26A...697A.106B/exportcitation).

# Contact
Email: Alexandre.Bougakov@obspm.fr
