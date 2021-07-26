function r8mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% r8mat_write writes an R8MAT file.
%
%  Discussion:
%
%    An R8MAT is an array of R8's.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    string OUTPUT_FILENAME, the output filename.
%
%    integer M, the spatial dimension.
%
%    integer N, the number of points.
%
%    real TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'R8MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
%  Alternative print statements:
%
%     fprintf ( output_unit, '  %14.6e', table(i,j) );
%     fprintf ( output_unit, '  %24.16e', table(i,j) );
%
  for i = 1 : m
    for j = 1 : n
      fprintf ( output_unit, '  %g', table(i,j) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
