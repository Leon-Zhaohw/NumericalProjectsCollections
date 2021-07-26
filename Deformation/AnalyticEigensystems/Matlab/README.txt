These Matlab scripts numerically verify the of the eigenpairs we derived.

The expressions that are output may appear slightly different from those listed
in the paper because they were hand-simplified to be more readable, and
constants such as mu and kappa were pinned to one to provide additional aid
to the simplification routine.

In order for these scripts to run, you MUST have either:
 - Matlab with the Symbolic Math Toolbox
 - Octave with the 'Symbolic' package (https://octave.sourceforge.io/symbolic/)

Verification scripts are provided for four different energies:
 - ARAP
 - Corotational
 - Symmetric ARAP
 - Symmetric Dirichlet

For each test, the PK1 is verified by performing a finite difference
convergence test against the original energy, and the Hessian is verified by
then performing an equivalent test against the PK1.

A symbolic version of the Hessian is then constructed, as well as symbolic
versions of the eigenpairs. The eigenpairs are build using the generic
formulas presented in the paper.

The eigenvector for each pair is then multipled against the symbolic Hessian, 
and it is verified that the result is the same as the proposed eigenvector 
scaled by the proposed eigenvalue.

The 2D and 3D cases are in separate subdirectories. To run all the tests in
2D, execute:
  
  cd 2D
  Verify_2D

To run the 3D tests, execute:

  cd 3D
  Verify_3D

If you try to run one after the other, you may get an error because some 2D
and 3D function names are duplicated. If this occurs, restart Octave/Matlab so 
that the search path is flushed and try again.
