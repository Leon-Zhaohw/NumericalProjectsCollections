function dg1d_advection_test ( )

%*****************************************************************************80
%
%% DG1D_ADVECTION_TEST tests DG1D_ADVECTION.
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
%    04 January 2019
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
  addpath ( '../dg1d_advection' )

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, ' DG1D_ADVECTION_TEST:\n' );
  fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Test DG1D_ADVECTION.\n' );

  Globals1D;
%
%  Order of polymomials used for approximation.
%
  N = 8;
  Elements = 10;
%
%  Generate simple mesh.
%
  [ Nv, VX, K, EToV ] = MeshGen1D ( 0.0, 2.0, Elements );
%
%  Initialize solver and construct grid and metric
%
  StartUp1D ( );
%
%  Set initial conditions.
%
  u = sin ( x );
%
%  Solve the Problem.
%
  FinalTime = 10.0;
  u = Advec1D ( u, FinalTime );

  filename = 'dg1d_advection_test.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, ' DG1D_ADVECTION:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../dg1d_advection' )

  return
end

