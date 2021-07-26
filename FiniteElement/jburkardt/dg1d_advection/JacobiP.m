function [ P ] = JacobiP ( x, alpha, beta, N )

%*****************************************************************************80
%
%% JACOBIP evaluates Jacobi polynomials.
%
%  Discussion:
%
%    (alpha,beta) > -1 and (alpha+beta <> -1).
%
%    The polynomials are evaluated at points x for order N and returned
%    in P[1:length(xp))]
%
%    The polynomials are normalized to be orthonormal.
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

%
%  Turn points into row if needed.
%
  xp = x; 
  dims = size(xp);
  if (dims(2)==1) 
    xp = xp'; 
  end;

  PL = zeros(N+1,length(xp)); 
%
%  Initial values P_0(x) and P_1(x)
%
  gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
    gamma(beta+1)/gamma(alpha+beta+1);
  PL(1,:) = 1.0/sqrt(gamma0);

  if (N==0) 
    P=PL';
    return;
  end;

  gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
  PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);

  if (N==1) 
    P=PL(N+1,:)'; 
    return;
  end;
%
%  Repeat value in recurrence.
%
  aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));
%
%  Forward recurrence using the symmetry of the recurrence.
%
  for i=1:N-1
    h1 = 2*i+alpha+beta;
    anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*...
      (i+1+beta)/(h1+1)/(h1+3));
    bnew = - (alpha^2-beta^2)/h1/(h1+2);
    PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
    aold =anew;
  end;

  P = PL(N+1,:)';

  return
end
