function dg1d_poisson_test ( )

%*****************************************************************************80
%
%% DG1D_POISSON_TEST tests the DG1D_POISSON library.
%
%  Discussion:
%
%   The formula for the value of the computed solution at the point
%   x in subinterval i, using quadratic polyonomials, is:
%
%     uh(x) = c[1+(i-1)*locdim] * dg1d_poisson_monomial(x,i,ne,0) 
%           + c[2+(i-1)*locdim] * dg1d_poisson_monomial(x,i,ne,1) 
%           + c[3+(i-1)*locdim] * dg1d_poisson_monomial(x,i,ne,2)
%
%  Modified:
%
%    06 January 2019
%
%  Author:
%
%    John Burkardt
%
  addpath ( '../dg1d_poisson' );

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_POISSON_TEST:\n' );
  fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Test dg1d_poisson\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  DG1D_POISSON applies the Discontinuous Galerkin Method (DG)\n' );
  fprintf ( 1, '  to a 1D Poisson problem.\n' );

  nel = 4;
  h = 1.0 / nel;
  ss = 1.0;
  penal = 1.0;
%
%  LOCDIM is actually the number of monomials used in each subinterval.
%
  locdim = 3;

  fprintf ( 1, '\n' );
  fprintf ( 1, '  0.0 < x < 1.0\n' );
  fprintf ( 1, '  Number of subintervals = %d\n', nel );
  fprintf ( 1, '  Number of monomials in expansion = %d\n', locdim );
  fprintf ( 1, '  Penalty parameter = %g\n', penal );
  fprintf ( 1, '  DG choice = %g\n', ss );
%
%  Compute C, which represents the solution as a set of coefficients
%  of monomials, indexed by polynomial degree and by subinterval.
%
  c = dg1d_poisson ( nel, ss, penal, @dg1d_poisson_test_source );
%
%  Evaluate the computed solution at M=5 points in each subinterval.
%
  m = 5;

  xh = zeros(m*nel,1);
  uh = zeros(m*nel,1);

  order = locdim;
  k = 0;
  for i = 1 : nel
    xl = ( i - 1 ) / nel;
    xr =   i       / nel;
    for j = 0 : m - 1
      k = k + 1;
      xh(k) = ( ( m - 1 - j ) * xl + j * xr ) / ( m - 1 );
      uh(k) = dg1d_poisson_interp ( xh(k), i, nel, order, c );
    end
  end
%
%  Tabulate the exact and computed solutions.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, '  I      X(I)     U(X(I))    Uh(X(I))\n' );
  fprintf ( 1, '\n' );
  for k = 1 : m * nel
    exact = dg1d_poisson_test_exact ( xh(k) );
    fprintf ( 1, '%2d  %8f  %8f  %8f\n', k, xh(k), exact, uh(k) );
  end
%
%  Evaluate the true solution at lots of points.
%
  x = linspace ( 0.0, 1.0, 101 );
  u = dg1d_poisson_test_exact ( x );
%
%  Make a plot comparing the exact and computed solutions.
%
  clf
  hold on
  plot ( xh, uh )
  plot ( x, u )
  grid ( 'on' )
  xlabel ( '<---X--->' );
  ylabel ( '<---U(X)--->' );
  title ( 'Compare exact and approximate solutions' )
  legend ( 'approximate', 'exact' )
  hold off
%
%  Save the plot to a file.
%
  filename = 'dg1d_poisson_test.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Graphics information saved in "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_POISSON_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../dg1d_poisson' );

  return
end

