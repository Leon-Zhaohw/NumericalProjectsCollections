function value = dg1d_poisson_monomial ( x, i, nel, order )

%*****************************************************************************80
%
%% DG1D_POISSON_MONOMIAL evaluates a monomial for the Poisson problem.
%
%  Discussion:
%
%    Assume the I-th interval has width H and midpoint XM.
%    Then the monomial of order ORDER is:
%
%      m(x;i,order) = ( 2 * ( x - xm ) / h ) ^ order
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
%    integer I, the index of the interval.
%
%    integer NEL, the number of intervals.
%
%    integer ORDER, the order of the monomial.
%
%  Output:
%
%    real VALUE, the value of the monomial at X.
%
  h = 1.0 / nel;
  xl = ( i - 1 ) * h;
  xr =   i       * h;
  xm = 0.5 * ( xl + xr );
  value = ( 2.0 * ( x - xm ) / h ) .^ order;

  return
end

