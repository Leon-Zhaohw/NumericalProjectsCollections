function [ rhsu ] = HeatCstabRHS1D ( u, time )

%*****************************************************************************80
%
%% HEATCSTABRHS1D: RHS flux in 1D heat equation using stabilized central flux.
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
  Globals1D;
%
%  Define field differences at faces
%
  du   = zeros(Nfp*Nfaces,K);
  du(:) = (u(vmapM)-u(vmapP))/2.0;
%
%  Impose boundary condition -- Dirichlet BC's
%
  uin  = -u(vmapI);
  du (mapI) = (u(vmapI) -  uin )/2.0;
  uout = -u(vmapO);
  du (mapO) = (u(vmapO) - uout)/2.0;
%
%  Compute q
%
  q = rx.*(Dr*u) - LIFT*(Fscale.*(nx.*du));
  dq = zeros(Nfp*Nfaces,K);
  dq(:) = q(vmapM)-q(vmapP);
%
%  Impose boundary condition -- Neumann BC's
%
  qin = q(vmapI);
  dq (mapI) = q(vmapI)-  qin; 
  qout = q(vmapO);
  dq (mapO) = q(vmapO)- qout;
%
%  Evaluate fluxes
%
  tau = 1.0;
  fluxq = nx.*(dq/2.0+tau*nx.*du);
%
%  Compute right hand sides of the semi-discrete PDE
%
  rhsu = rx.*(Dr*q) - LIFT*(Fscale.*fluxq);

  return
end

