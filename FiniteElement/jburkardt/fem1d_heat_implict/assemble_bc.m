function [ a, b ] = assemble_bc ( x_num, x, t, bc_fun, a, b )

%*****************************************************************************80
%
%% ASSEMBLE_BC modifies the linear system for Dirichlet boundary conditions.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 February 2012
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer X_NUM, the number of nodes.
%
%    real X(X_NUM), the coordinates of nodes.
%
%    real T, the current time.
%
%    real BC_FUN(), a function to set the Dirichlet conditions.
%
%    sparse real A(X_NUM,X_NUM), the coefficient matrix.
%
%    real B(X_NUM), the right hand side.
%
%  Output:
%
%    sparse real A(X_NUM,X_NUM), the adjusted matrix.
%
%    real B(X_NUM), the adjusted right hand side.
%
  u = zeros ( x_num, 1 );
  u(1:x_num,1) = bc_fun ( x_num, x, t, u );

  a(1,1:x_num) = 0.0;
  a(1,1) = 1.0;
  b(1) = u(1);

  a(x_num,1:x_num) = 0.0;
  a(x_num,x_num) = 1.0;
  b(x_num) = u(x_num);

  return
end
