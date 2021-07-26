function value = dg1d_poisson_test_source ( x )

%*****************************************************************************80
%
%% DG1D_POISSON_TEST_SOURCE evaluates the source function for the test example.
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
%    real VALUE, the value of the source term at X.
%
  value = - ( 2.0 * x - 2.0 * ( 1.0 - 2.0 * x ) ...
    + 4.0 * x * ( x - x * x ) ) * exp ( - x * x );

  return
end

