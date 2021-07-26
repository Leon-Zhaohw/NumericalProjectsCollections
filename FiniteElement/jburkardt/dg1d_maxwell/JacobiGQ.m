function [ x, w ] = JacobiGQ ( alpha, beta, N )

%*****************************************************************************80
%
%% JACOBIGQ computes Nth order Gauss quadrature points and weights.
%
%  Discussion:
%
%    The N'th order Gauss quadrature points, x, and weights, w, are associated 
%    with the Jacobi polynomial, of type (alpha,beta) > -1 ( <> -0.5).
%
%  Licensing:
%
%    Permission to use this software for noncommercial
%    research and educational purposes is hereby granted
%    without fee.  Redistribution, sale, or incorporation
%    of this software into a commercial product is prohibited.
%
%    THE AUTHORS OR PUBLISHER DISCLAIMS ANY AND ALL WARRANTIES
%    WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
%    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
%    PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHORS OR
%    THE PUBLISHER BE LIABLE FOR ANY SPECIAL, INDIRECT OR
%    CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
%    RESULTING FROM LOSS OF USE, DATA OR PROFITS.
%
%  Modified:
%
%    17 September 2018
%
%  Author:
%
%    Original version by Jan Hesthaven, Tim Warburton.
%    Some modifications by John Burkardt.
%
%  Reference:
%
%    Jan Hesthaven, Tim Warburton,
%    Nodal Discontinuous Galerkin Methods: 
%    Algorithms, Analysis, and Applications,
%    Springer, 2007,
%    ISBN: 978-0387720654.
%
  if (N==0) 
    x(1)= -(alpha-beta)/(alpha+beta+2); 
    w(1) = 2.0; 
    return;
  end;
%
%  Form symmetric matrix from recurrence.
%
  J = zeros(N+1);
  h1 = 2*(0:N)+alpha+beta;
  J = diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1) + ...
    diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta).*...
    ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);

  if (alpha+beta<10*eps) 
    J(1,1) = 0.0;
  end;

  J = J + J';
%
%  Compute quadrature by eigenvalue solve
%
  [V,D] = eig(J); 
  x = diag(D);

  w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
    gamma(beta+1)/gamma(alpha+beta+1);

  return;
end
