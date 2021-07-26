function fem1d_bvp_linear_test00 ( )

%*****************************************************************************80
%
%% FEM1D_BVP_LINEAR_TEST00 tests FEM1D_BVP_LINEAR.
%
%  Discussion:
%
%    - uxx + u = x  for 0 < x < 1
%    u(0) = u(1) = 0
%
%    exact  = x - sinh(x) / sinh(1)
%    exact' = 1 - cosh(x) / sinh(1)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Dianne O'Leary,
%    Scientific Computing with Case Studies,
%    SIAM, 2008,
%    ISBN13: 978-0-898716-66-5,
%    LC: QA401.O44.
%
  n = 11;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'FEM1D_BVP_LINEAR_TEST00\n' );
  fprintf ( 1, '  Solve -( A(x) U''(x) )'' + C(x) U(x) = F(x)\n' );
  fprintf ( 1, '  for 0 < x < 1, with U(0) = U(1) = 0.\n' );
  fprintf ( 1, '  A(X)  = 1.0\n' );
  fprintf ( 1, '  C(X)  = 1.0\n' );
  fprintf ( 1, '  F(X)  = X\n' );
  fprintf ( 1, '  U(X)  = X - SINH(X) / SINH(1)\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of nodes = %d\n', n );
%
%  Geometry definitions.
%
  x_first = 0.0;
  x_last = 1.0;
  x = linspace ( x_first, x_last, n );
  x = x(:);

  u = fem1d_bvp_linear ( n, @a00, @c00, @f00, x );

  uexact = exact00 ( x );

  fprintf ( 1, '\n' );
  fprintf ( 1, '     I    X         U         Uexact    Error\n' );
  fprintf ( 1, '\n' );

  for i = 1 : n
    fprintf ( 1, '  %4d  %8f  %8f  %8f  %8e\n', ...
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) ) );
  end

  e1 = l1_error ( n, x, u, @exact00 );
  e2 = l2_error_linear ( n, x, u, @exact00 );
  h1s = h1s_error_linear ( n, x, u, @exact_ux00 );
  mx = max_error_linear ( n, x, u, @exact00 );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  l1 norm of error  = %g\n', e1 );
  fprintf ( 1, '  L2 norm of error  = %g\n', e2 );
  fprintf ( 1, '  Seminorm of error = %g\n', h1s );
  fprintf ( 1, '  Max norm of error = %g\n', mx );

  return
end
function value = a00 ( x )

%*****************************************************************************80
%
%% A00 evaluates A function #0.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real X, the evaluation point.
%
%  Output:
%
%    real VALUE, the value of A(X).
%
  value = 1.0;

  return
end
function value = c00 ( x )

%*****************************************************************************80
%
%% C00 evaluates C function #0.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real X, the evaluation point.
%
%  Output:
%
%    real VALUE, the value of C(X).
%
  value = 1.0;

  return
end
function value = exact00 ( x )

%*****************************************************************************80
%
%% EXACT00 evaluates exact solution #0.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real X, the evaluation point.
%
%  Output:
%
%    real VALUE, the value of U(X).
%
  value = x - sinh ( x ) / sinh ( 1.0 );

  return
end
function value = exact_ux00 ( x )

%*****************************************************************************80
%
%% EXACT_UX00 evaluates the derivative of exact solution #0.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real X, the evaluation point.
%
%  Output:
%
%    real VALUE, the value of dUdX(X).
%
  value = 1.0 - cosh ( x ) / sinh ( 1.0 );

  return
end
function value = f00 ( x )

%*****************************************************************************80
%
%% F00 evaluates right hand side function #0.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real X, the evaluation point.
%
%  Output:
%
%    real VALUE, the value of F(X).
%
  value = x;

  return
end

