function [ x ] = JacobiGL ( alpha, beta, N )

%*****************************************************************************80
%
%% JACOBIGL computes the N'th order Gauss Lobatto quadrature points.
%
%  Discussion:
%
%    These points are associated with the Jacobi polynomial,
%    of type (alpha,beta) > -1 ( <> -0.5). 
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
  x = zeros(N+1,1);

  if ( N == 1 )
    x(1)=-1.0;
    x(2)=1.0; 
    return; 
  end;

  [xint,w] = JacobiGQ(alpha+1,beta+1,N-2);

  x = [ -1, xint', 1 ]';

  return;
end
