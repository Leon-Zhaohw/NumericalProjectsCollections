
This Matlab package is developed by Corinna Bahriawati
(corinna@aurora.anum.tuwien.ac.at) and Carsten Carstensen (cc@math.hu-berlin.de)
dated Feb, 2005. This package will be released through ACM Transaction
on Mathematical Software. The author is not responsible for any damage
that may be related to using this package.

1. General information on the package

There are three directories when the package is unzipped:

      Documentation
      CBaCCMFEM
      Examples

There are two papers in "documentation" directory:

  [1]. "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
        MFEM with a Posteriori Error Control", Corinna Bahriawati and 
        Carsten Carstensen, submitted to ACM Trans. Math. Software.

  [2]. "User Manual Three Matlab Implementations of the Lowest-Order Raviart-Thomas
        MFEM with a Posteriori Error Control", Corinna Bahriawati and 
        Carsten Carstensen, submitted to ACM Trans. Math. Software.

The paper [1] presents the overall theory, while [2]
describes the software package and the examples in detail.

The directory "CBaCCMFEM" contains the Matlab modules of "Three Matlab
Implementations of the Lowest-Order Raviart-Thomas MFEM with a Posteriori Error Control
package".

The directory "examples" contains the Matlab modules and data that generate
numerical examples. [The method and application
only treat 2-dimensional problem.]

2. Installation

Unzip and save all the files, open Matlab, set proper path to
the directory containing the Matlab modules, the package is
now ready to run.

3. Simple applications

To compute the numerical approximation of the Poisson equation with
inhomogeneous mixed boundary conditions in 2D with lowest-order Raviart-Thomas
mixed finite elements in Matlab format, simply call

>> EBmfem

Here we employ the edge-basis functions for the lowest-order
Raviart Thomas finite element.
The output contains the vector of the numerical solution of the displacement $x$ resp. the flux $p$. 
The first element number of the vector list discrete displacements,
and the  shows discrete fluxes.


4. Applications

There are three main programs:

1. EBmfem:    calculates the numerical solution of the problem
              based on the edge-basis functions for the lowest-order 
              Raviart Thomas finite element.


2. LMmfem:   calculates the numerical solution of the problem
              based on the Lagrange multiplier technique for the lowest-order 
              Raviart Thomas finite element.

       
3. CRmfem:   calculates the numerical solution of the problem
              based on the non-conforming Crouzeix-Raviart for the lowest-order 
              Raviart Thomas finite element.

Users should read the papers attached in the package. At least
the paper [2] above. All modules contain "help" message that can be
accessed in Matlab. You may try

>> help EBmfem



