This Matlab script symbolically verifies the eigenpairs we have derived for IV_C, a.k.a. I_5.

In order for these scripts to run, you MUST have either:
 - Matlab with the Symbolic Math Toolbox
 - Octave with the 'Symbolic' package (https://octave.sourceforge.io/symbolic/)

From the top directory, invoke:

    Verify_I5

The script will numerically verify our PK1 and Hessian expressions, and then symbolically verify our eigenpairs. It will also verify that the system is rank three, so there are no lingering eigenpairs that we missed.