function value = dg1d_poisson_test_exact ( x )

%*****************************************************************************80
%
%% DG1D_POISSON_TEST_EXACT evaluates the exact solution.
%
%  Modified:
%
%    15 September 2018
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real X, the argument.
%
%  Output:
%
%    real VALUE, the exact solution at X.
%
  value = ( 1.0 - x ) .* exp ( - x .^ 2 );

  return
end

