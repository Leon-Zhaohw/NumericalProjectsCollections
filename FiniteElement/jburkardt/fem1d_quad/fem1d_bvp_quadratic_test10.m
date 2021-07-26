function fem1d_bvp_quadratic_test10 ( )

%*****************************************************************************80
%
%% FEM1D_BVP_QUADRATIC_TEST10 tests FEM1D_BVP_LINEAR.
%
%  Discussion:
%
%    We want to compute errors and do convergence rates for the 
%    following problem:
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
%    19 July 2015
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
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FEM1D_BVP_QUADRATIC_TEST10\n' );
  fprintf ( 1, '  Solve -( A(x) U''(x) )'' + C(x) U(x) = F(x)\n' );
  fprintf ( 1, '  for 0 < x < 1, with U(0) = U(1) = 0.\n' );
  fprintf ( 1, '  A(X)  = 1.0\n' );
  fprintf ( 1, '  C(X)  = 1.0\n' );
  fprintf ( 1, '  F(X)  = X\n' );
  fprintf ( 1, '  U(X)  = X - SINH(X) / SINH(1)\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, ' log(E)    E         L2error         H1error        Maxerror\n' );
  fprintf ( 1, '\n' );

  e_log_max = 6;

  ne_plot = zeros ( e_log_max + 1 );
  h_plot = zeros ( e_log_max + 1 );
  l2_plot = zeros ( e_log_max + 1 );
  h1_plot = zeros ( e_log_max + 1 );
  mx_plot = zeros ( e_log_max + 1 );

  for e_log = 0 : e_log_max

    ne = 2 ^ ( e_log + 1 );
    n = ne + 1;

    x_lo = 0.0;
    x_hi = 1.0;
    x = linspace ( x_lo, x_hi, n );

    u = fem1d_bvp_quadratic ( n, @a10, @c10, @f10, x );

    ne_plot(e_log+1) = ne;

    h_plot(e_log+1) = ( x_hi - x_lo ) / ne;

    l2_plot(e_log+1) = l2_error_quadratic ( n, x, u, @exact10 );

    h1_plot(e_log+1) = h1s_error_quadratic ( n, x, u, @exact_ux10 );

    mx_plot(e_log+1) = max_error_quadratic ( n, x, u, @exact10 );

    fprintf ( 1, '  %4d  %4d  %14g  %14g  %14g\n', ...
      e_log, ne, l2_plot(e_log+1), h1_plot(e_log+1), mx_plot(e_log+1) );

  end

  fprintf ( 1, '\n' );
  fprintf ( 1, ' log(E1)  E1 / E2          L2rate          H1rate         Maxrate\n' );
  fprintf ( 1, '\n' );

  for e_log = 0 : e_log_max - 1
    ne1 = ne_plot(e_log+1);
    ne2 = ne_plot(e_log+2);
    r = ne2 / ne1;
    l2 = l2_plot(e_log+1) / l2_plot(e_log+2);
    l2 = log ( l2 ) / log ( r );
    h1 = h1_plot(e_log+1) / h1_plot(e_log+2);
    h1 = log ( h1 ) / log ( r );
    mx = mx_plot(e_log+1) / mx_plot(e_log+2);
    mx = log ( mx ) / log ( r );
    fprintf ( 1, '  %4d  %4d/%4d  %14g  %14g  %14g\n', ...
      e_log, ne1, ne2, l2, h1, mx );
  end
%
%  Plot the L2 error as a function of NE.
%
  figure ( );
  clf
  loglog ( ne_plot, l2_plot, 'bo-' );
  xlabel ( '<---NE--->' )
  ylabel ( '<---L2(error)--->' )
  title ( 'L2 error as function of number of elements' )
  grid on
  axis 'equal'
  filename = 'l2error_test10.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Created plot file "%s".\n', filename );
%
%  Plot the max error as a function of NE.
%
  figure ( );
  loglog ( ne_plot, mx_plot, 'bo-' );
  xlabel ( '<---NE--->' )
  ylabel ( '<---Max(error)--->' )
  title ( 'Max error as function of number of elements' )
  grid on
  axis equal
  filename = 'maxerror_test10.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Created plot file "%s".\n', filename );
%
%  Plot the H1 error as a function of NE.
%  Once again MATLAB refuses to do what I say.
%  AXIS EQUAL, AXIS EQUAL, AXIS EQUAL!
%
  figure ( );
  loglog ( ne_plot, h1_plot, 'bo-' );
  axis equal
  xlabel ( '<---NE--->' )
  ylabel ( '<---H1(error)--->' )
  title ( 'H1 error as function of number of elements' )
  grid on

  filename = 'h1error_test10.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Created plot file "%s".\n', filename );

  return
end
function value = a10 ( x )

%*****************************************************************************80
%
%% A10 evaluates A function #10.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 July 2015
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
function value = c10 ( x )

%*****************************************************************************80
%
%% C10 evaluates C function #10.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 July 2015
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
function value = exact10 ( x )

%*****************************************************************************80
%
%% EXACT10 evaluates exact solution #10.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 July 2015
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
function value = exact_ux10 ( x )

%*****************************************************************************80
%
%% EXACT_UX10 evaluates the derivative of exact solution #10.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 July 2015
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
function value = f10 ( x )

%*****************************************************************************80
%
%% F10 evaluates right hand side function #10.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 July 2015
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

