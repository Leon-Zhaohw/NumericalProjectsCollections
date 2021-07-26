function dg1d_heat_test ( )

%*****************************************************************************80
%
%% DG1D_HEAT_TEST tests DG1D_HEAT.
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
%    06 January 2019
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
  addpath ( '../dg1d_heat' );

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_HEAT_TEST:\n' );
  fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Test dg1d_heat\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Discontinuous Galerkin method for the 1D heat equation.\n' );
%
%  Run Globals1D script to define global variables.
%
  Globals1D;
%
%  Polynomial order used for approximation.
%
  N = 8;
%
%  Set up the Mesh.
%
  [ Nv, VX, K, EToV ] = MeshGen1D ( 0, 2*pi, 20 );
%
%  Run StartUp1D script to initialize solver and construct grid and metric.
%
  StartUp1D;
%
%  Set the initial conditions.
%
  u = sin ( x );
%
%  Solve the problem.
%
  FinalTime = 0.8;
  [ u, time ] = Heat1D ( u, FinalTime );
%
%  Rearrange the data as vectors.
%
  xv = reshape ( x, (N+1)*20, 1 );
  uv = reshape ( u, (N+1)*20, 1 );
%
%  Tabulate the solution.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, '     I      X(I)      U(I)\n' );
  fprintf ( 1, '\n' );
  for k = 1 : ( N + 1 ) * 20
    fprintf ( 1, '  %4d  %8f  %8f\n', k, xv(k), uv(k) );
  end
%
%  Plot the solution.
%
  plot ( xv, uv );
  grid ( 'on' );
  xlabel ( '<---X--->' )
  ylabel ( '<--U(X,T)-->' )
  title ( 'Heat Equation Solution at Final Time' )
%
%  Save the plot in a file.
%
  filename = 'dg1d_heat_test.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Graphics saved in file "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_HEAT_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../dg1d_heat' );

  return
end
