function dg1d_maxwell_test ( )

%*****************************************************************************80
%
%% DG1D_MAXWELL_TEST tests DG1D_MAXWELL.
%
%  Licensing:
%
%    Permission to use this software for noncommercial
%    research and educational purposes is hereby granted
%    without fee.  Redistribution, sale, or incorporation
%    of this software into a commercial product is prohibited.
%
%    THE AUTHORS OR PUBLISHER DISCLAIMS ANY AND ALL WARRANTIES
%    WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
%    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
%    PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHORS OR
%    THE PUBLISHER BE LIABLE FOR ANY SPECIAL, INDIRECT OR
%    CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
%    RESULTING FROM LOSS OF USE, DATA OR PROFITS.
%
%  Modified:
%
%    24 September 2018
%
%  Author:
%
%    Original version by Jan Hesthaven, Tim Warburton.
%    Some modifications by John Burkardt.
%
%  Reference:
%
%    Jan Hesthaven, Tim Warburton,
%    Nodal Discontinuous Galerkin Methods: 
%    Algorithms, Analysis, and Applications,
%    Springer, 2007,
%    ISBN: 978-0387720654.
%
  addpath ( '../dg1d_maxwell' );

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_MAXWELL_TEST:\n' );
  fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Test dg1d_maxwell\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Solve the 1D Maxwell equations using the discontinuous\n' );
  fprintf ( 1, '  Galerkin method.\n' );
  
  Globals1D;
%
%  Polynomial order used for approximation.
%
  N = 6;
%
%  Generate simple mesh.
%
  [ Nv, VX, K, EToV ] = MeshGen1D ( -2.0, 2.0, 80 );
%
%  Initialize solver and construct grid and metric.
%
  StartUp1D;
%
%  Set up material parameters.
%
  eps1 = [ ones(1,K/2), 2*ones(1,K/2) ];
  mu1 = ones(1,K);
  epsilon = ones(Np,1)*eps1;
  mu = ones(Np,1)*mu1;
%
%  Set initial conditions.
%
  E = sin(pi*x) .* (x<0);
  H = zeros(Np,K);
%
%  Solve Problem.
%
  FinalTime = 10;
  [ E, H ] = Maxwell1D ( E, H, epsilon, mu, FinalTime );
%
%  Rearrange the data as vectors.
%
  xv = reshape ( x, (N+1)*80, 1 );
  ev = reshape ( E, (N+1)*80, 1 );
  hv = reshape ( H, (N+1)*80, 1 );
%
%  Tabulate the solution.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, '     I      X(I)      E(I)          H(I)\n' );
  fprintf ( 1, '\n' );
  for k = 1 : ( N + 1 ) * 80
    fprintf ( 1, '  %4d  %8f  %8f  %8f\n', k, xv(k), ev(k), hv(k) );
  end
%
%  Plot the solutions.
%
  figure ( 1 )
  plot ( xv, ev, 'r-', xv, hv, 'b-' );
  grid ( 'on' );
  xlabel ( '<---X--->' )
  ylabel ( '<--E(X,T)(red), H(X,T)(blue)-->' )
  title ( 'Maxwell Equation Solution at Final Time' )
%
%  Save the plot in a file.
%
  filename = 'dg1d_maxwell_test.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Graphics saved in file "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_MAXWELL_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  return
end

 
