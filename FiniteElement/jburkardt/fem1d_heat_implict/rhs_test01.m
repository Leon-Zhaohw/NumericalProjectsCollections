function rhs_value = rhs_test01 ( x_num, x, t )

%*****************************************************************************80
%
%% RHS_TEST01 evaluates the right hand side function for problem 1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    29 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer X_NUM, the number of evaluation points.
%
%    real X(X_NUM,1), the evaluation points.
%
%    real T, the time.
%
%  Output:
%
%    real RHS_VALUE(X_NUM,1), the right hand side function at
%    the given positions and time T.
%
  rhs_value(1:x_num,1) = zeros ( x_num, 1 );

  return
end
