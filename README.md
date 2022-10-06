# stokes-dg-experiments
Script that runs experiments on the implemented H(div)-conforming elements for Stokes equations in ParMooN.

## Experiments
Simulations with the ParMooN package to numerically solve the Stokes equations using the finite element method with various discretisation types
and function spaces. It was used for studying the properties of the Discontinous Galerkin (DG) implementation of the Stokes equations in ParMooN,
to which I contributed as part of my work at WIAS Berlin and my master thesis at FU Berlin.

 Two steady state examples are used to evaluate the performance of the numerical solvers, with analytic solutions for pairs of velocity and pressure 
 defined as "harmonic" and "polynomial". (For full definitions and visualisations see pages 38 - 40 in the [master thesis](https://github.com/cristina-v-melnic/stokes-dg-experiments/blob/main/Master_Thesis_signed.pdf).)

### Code functionality
On these examples, three kinds of numerical experiments can be performed by running the functions below in main().
- *loop_refinement_all_examples()* - Generate a set of convergence histories with different examples, velocity spaces and discretization
types in order to compare
  the optimality of different methods.
- *get_Re_nr_scaling()* - Function to study pressure robustness by generating sets of convergence histories with different Reynolds number 
and observing whether or not the velocity error scales with $Re = \frac{1}{\nu}$. 
- *check_optimal_sigma()* - Determine the optimal choice of the sigma parameter for the DG method by comparing the magnitudes of errors 
from simulations across ranges of values for sigma (grid-search). 

The analysis of the generated datasets by the simulations above is contained in the plots made with the [stokes-dg-figures project](https://github.com/cristina-v-melnic/stokes-dg-figures). 
(For discussions and graphics see sections 4.2.- 4.4. and the appendix of the [master thesis](https://github.com/cristina-v-melnic/stokes-dg-experiments/blob/main/Master_Thesis_signed.pdf).)

## Credits
This work was part of my master's thesis project with Prof. Dr. Volker John, Derk Frerichs and Dr. Ulrich Wilbrandt from the
"Numerical Mathematics and Scientific Computing" group at the Weierstrass Institute Berlin. It was financially supported by
the DAAD funding of my master's studies, and the "WIAS Female Master Students Program" funding of my research assistanship.
