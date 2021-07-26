function fem1d_spectral_numeric_test ( )

%*****************************************************************************80
%
%% FEM1D_SPECTRAL_NUMERICY_TEST tests FEM1D_SPECTRAL_NUMERIC.
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
  addpath ( '../fem1d_spectral_numeric' );

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FEM1D_SPECTRAL_NUMERIC_TEST:\n' );
  fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Test FEM1D_SPECTRAL_NUMERIC.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Run FEM1D_SPECTRAL_NUMERIC with an increasing number of basis functions.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '   N        L2-error        H1-error\n' );
  fprintf ( 1, '\n' );

  for n = 1 : 2 : 11
    [ l2, h1 ] = fem1d_spectral_numeric ( n );
    fprintf ( 1, '  %2d  %14.6g  %14.6g\n', n, l2, h1 );
  end
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FEM1D_SPECTRAL_NUMERIC_TEST:\n' );
  fprintf ( 1, '  Normal end of execution\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../fem1d_spectral_numeric' );

  return
end

