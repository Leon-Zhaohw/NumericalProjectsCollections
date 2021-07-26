function fem1d_nonlinear_test ( )

%*****************************************************************************80
%
%% FEM1D_NONLINEAR_TEST tests FEM1D_NONLINEAR.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 January 2019
%
%  Author:
%
%    John Burkardt
%
  addpath ( '../fem1d_nonlinear' )

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FEM1D_NONLINEAR_TEST:\n' );
  fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Test FEM1D_NONLINEAR.\n' );

  problem = 1;
  fem1d_nonlinear ( problem );

  problem = 2;
  fem1d_nonlinear ( problem );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FEM1D_NONLINEAR_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../fem1d_nonlinear' )

  return
end
