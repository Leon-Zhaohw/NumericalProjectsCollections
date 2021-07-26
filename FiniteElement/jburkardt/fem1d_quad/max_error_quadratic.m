function value = max_error_quadratic ( n, x, u, exact, value )

%*****************************************************************************80
%
%% MAX_ERROR_QUADRATIC estimates the max error norm of a finite element solution.
%
%  Discussion:
%
%    We assume the finite element method has been used, over an interval [A,B]
%    involving N nodes, with piecewise quadratic elements used for the basis.
%    The coefficients U(1:N) have been computed, and a formula for the
%    exact solution is known.
%
%    This function estimates the max norm of the error:
%
%      MAX_NORM = Integral ( A <= X <= B ) max ( abs ( U(X) - EXACT(X) ) ) dX
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer N, the number of nodes.
%
%    real X(N), the mesh points.
%
%    real U(N), the finite element coefficients.
%
%    function EQ = EXACT ( X ), returns the value of the exact
%    solution at the point X.
%
%  Output:
%
%    real VALUE, the estimated max norm of the error.
%
  quad_num = 8;
  value = 0.0;

  e_num = ( n - 1 ) / 2;

  for e = 1 : e_num

    l = 2 * e - 1;
    xl = x(l);
    ul = u(l);

    m = 2 * e;
    xm = x(m);
    um = u(m);

    r = 2 * e + 1;
    xr = x(r);
    ur = u(r);

    for q = 0 : quad_num - 1

      xq = ( ( quad_num - q ) * xl   ...
           + (            q ) * xr ) ...
         /   ( quad_num     );

      vl = ( ( xq - xm ) / ( xl - xm ) ) ...
         * ( ( xq - xr ) / ( xl - xr ) );

      vm = ( ( xq - xl ) / ( xm - xl ) ) ...
         * ( ( xq - xr ) / ( xm - xr ) );

      vr = ( ( xq - xl ) / ( xr - xl ) ) ...
         * ( ( xq - xm ) / ( xr - xm ) );

      uq = u(l) * vl + u(m) * vm + u(r) * vr;

      eq = exact ( xq );

      value = max ( value, abs ( uq - eq ) );

    end

  end
%
%  For completeness, check last node.
%
  xq = x(n);
  uq = u(n);
  eq = exact ( xq );

  value = max ( value, abs ( uq - eq ) );
%
%  Integral approximation requires multiplication by interval length.
%
  value = value * ( x(n) - x(1) );

  return
end
