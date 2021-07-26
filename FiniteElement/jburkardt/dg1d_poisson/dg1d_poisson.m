function c =  dg1d_poisson ( nel, ss, penal, f )

%*****************************************************************************80
%
%% DG1D_POISSON applies Discontinuous Galerkin method to 1D Poisson equation.
%
%  Discussion:
%
%    A 1D version of the Poisson equation has the form
%
%      - ( K(x) u'(x) )' = f(x)  for 0 < x < 1
%
%      u(0) = 1
%      u(1) = 0
%
%    Here, we will assume that K(x) = 1.
%
%    This function computes an approximate discrete solution to the problem,
%    using a version of the Discontinuous Galerkin method.  The interval
%    [0,1] is divided into NEL equal subintervals, over each of which a set
%    of LOCDIM=3 basis monomials are defined, centered at the midpoint of
%    the subinterval, and normalized to have unit value at the subinterval
%    endpoints.
%
%    The discontinous Galerkin equations are then set up as the linear
%    system A*c=b, where c represents coefficients of the basis functions.
%    The result of solving this system can then be used to tabulate, evaluate,
%    or plot the approximate DG solution function.
%
%  Modified:
%
%    17 September 2018
%
%  Author:
%
%    Original version by Beatrice Riviere.
%    Some modifications by John Burkardt.
%
%  Reference:
%
%    Beatrice Riviere,
%    Discontinuous Galerkin Methods for Solving Elliptic and Parabolic Equations,
%    SIAM, 2008
%    ISBN: 978-0-898716-56-6
%
%  Input:
%
%    integer NEL, the number of subintervals.
%
%    real SS, the symmetrization parameter.
%    Three values are meaningful:
%    1.0: NIPG, nonsymmetric interior penalty Galerkin method.
%    0.0: IIPG, incomplete interior penalty Galerkin method.
%   -1.0: SIPG, symmetric interior penalty Galerkin method.
%
%    real PENAL, the penalty parameter.
%
%    value = F(x), a function which evaluates the right
%    hand side of the Poisson equation.
%
%  Output:
%
%    real C(3*NEL), the DG coefficients.  The dimension of 3
%    is due to the use of piecewise quadratics in the interpolation scheme.
%

%
%  Dimension of local matrices.
%  This number should correspond to the number of monomials used in each
%  subinterval.  Because it is set to 3, we are using piecewise quadratics.
%
  locdim = 3;
%
%  Local matrices.
%  These matrices have order locdimxlocdim.
%  If we want to use higher order methods, then these local matrices would
%  need to be enlarged, and the appropriate additional values inserted.
%
  Amat = nel * [ ...
    0.0, 0.0, 0.0; ...
    0.0, 4.0, 0.0; ...
    0.0, 0.0, 16.0 / 3.0 ];

  Bmat = nel * [ ...
              penal,            1.0 - penal,            - 2.0 + penal; ...
       - ss - penal,     - 1.0 + ss - penal,         2.0 - ss - penal; ...
   2.0 * ss + penal, 1.0 - 2.0 * ss - penal, - 2.0 + 2.0 * ss + penal ];

  Cmat = nel * [ ...
               penal,            - 1.0 + penal,            - 2.0 + penal; ...
          ss + penal,       - 1.0 + ss + penal,       - 2.0 + ss + penal; ... 
    2.0 * ss + penal, - 1.0 + 2.0 * ss + penal, - 2.0 + 2.0 * ss + penal ];

  Dmat = nel * [ ...
               - penal,            - 1.0 + penal,            2.0 - penal; ...
          - ss - penal,       - 1.0 + ss + penal,       2.0 - ss - penal; ... 
    - 2.0 * ss - penal, - 1.0 + 2.0 * ss + penal, 2.0 - 2.0 * ss - penal ];

  Emat = nel * [ ...
               - penal,            1.0 - penal,            2.0 - penal; ...
            ss + penal,     - 1.0 + ss + penal,     - 2.0 + ss + penal; ...
    - 2.0 * ss - penal, 1.0 - 2.0 * ss - penal, 2.0 - 2.0 * ss - penal ];

  F0mat = nel * [ ...
                 penal,              2.0 - penal,            - 4.0 + penal; ... 
    - 2.0 * ss - penal, - 2.0 + 2.0 * ss + penal,   4.0 - 2.0 * ss - penal; ...
      4.0 * ss + penal,   2.0 - 4.0 * ss - penal, - 4.0 + 4.0 * ss + penal ];

  FNmat = nel * [ ...
               penal,            - 2.0 + penal,            - 4.0 + penal; ...
    2.0 * ss + penal, - 2.0 + 2.0 * ss + penal, - 4.0 + 2.0 * ss + penal; ...
    4.0 * ss + penal, - 2.0 + 4.0 * ss + penal, - 4.0 + 4.0 * ss + penal ];
%
%  Dimension of global matrix.
%
  glodim = nel * locdim;
%
%  Initialize the global data.
%
  A = zeros ( glodim, glodim );
  b = zeros ( glodim, 1 );
%
%  Gauss quadrature weights and points of order 2, which should be sufficient
%  for integrals of products of piecewise quadratic functions.
%
  ng = 2;
  wg = [ 1.0, 1.0 ];
  sg = [ -0.577350269189, 0.577350269189 ];
%
%  Assemble the global matrix and right hand side.
%

%
%  First subinterval.
%
  for ii = 1 : locdim
    for jj = 1 : locdim
      A(ii,jj) = A(ii,jj) + Amat(ii,jj) + F0mat(ii,jj) + Cmat(ii,jj);
      je = locdim + jj;
      A(ii,je) = A(ii,je) + Dmat(ii,jj);
    end
  end

  b(1) = nel * penal;
  b(2) = nel * penal * (-1.0) - ss * 2.0* nel;
  b(3) = nel * penal + ss * 4.0 * nel;

  for ig = 1 : ng
    xval = ( sg(ig) + 1.0 ) / ( 2.0 * nel );
    b(1) = b(1) + wg(ig) * f ( xval ) / ( 2.0 * nel ) * 1.0;
    b(2) = b(2) + wg(ig) * f ( xval ) / ( 2.0 * nel ) * sg(ig);
    b(3) = b(3) + wg(ig) * f ( xval ) / ( 2.0 * nel ) * sg(ig) * sg(ig);
  end
%
%  Intermediate subintervals.
%
  for i = 2 : nel - 1
    for ii = 1 : locdim
      ie = ii + ( i - 1 ) * locdim;
      for jj = 1 : locdim
        je = jj + ( i - 1 ) * locdim;
        A(ie,je) = A(ie,je) + Amat(ii,jj) + Bmat(ii,jj) + Cmat(ii,jj);
        je = jj + ( i - 2 ) * locdim;
        A(ie,je) = A(ie,je) + Emat(ii,jj);
        je = jj + i * locdim;
        A(ie,je) = A(ie,je) + Dmat(ii,jj);
      end
      for ig = 1 : ng
        xval = ( sg(ig) + 2.0 * ( i - 1 ) + 1.0 ) / ( 2.0 * nel );
        b(ie) = b(ie) + wg(ig) * f ( xval ) / ( 2.0 * nel ) * ( sg(ig) ^ ( ii - 1 ) );
      end
    end
  end
%
%  Last subinterval.
%
  for ii = 1 : locdim
    ie = ii + ( nel - 1 ) * locdim;
    for jj = 1 : locdim
      je = jj + ( nel - 1 ) * locdim;
      A(ie,je) = A(ie,je) + Amat(ii,jj) + FNmat(ii,jj) + Bmat(ii,jj);
      je = jj + ( nel - 2 ) * locdim;
      A(ie,je) = A(ie,je) + Emat(ii,jj);
    end
    for ig = 1 : ng
      xval = ( sg(ig) + 2.0 * ( nel - 1 ) + 1.0 ) / ( 2.0 * nel );
      b(ie) = b(ie) + wg(ig) * f ( xval ) / ( 2.0 * nel ) * ( sg(ig) ^ ( ii - 1 ) );
    end
  end
%
%  Solve the linear system.
%  
  c = A \ b;
    
  return
end

