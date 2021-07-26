function [ rhsE, rhsH ] = MaxwellRHS1D ( E, H, eps, mu )

%*****************************************************************************80
%
%% MAXWELLRHS1D evaluates the RHS flux in 1D Maxwell equations.
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
%    24 September 2018
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
%  Compute impedance.
%
  Zimp = sqrt ( mu ./ eps );
%
%  Define field differences at faces.
%
  dE = zeros(Nfp*Nfaces,K);
  dE(:) = E(vmapM)-E(vmapP);
  dH = zeros(Nfp*Nfaces,K);
  dH(:) = H(vmapM)-H(vmapP);
  Zimpm = zeros(Nfp*Nfaces,K);
  Zimpm(:) = Zimp(vmapM);
  Zimpp = zeros(Nfp*Nfaces,K);
  Zimpp(:) = Zimp(vmapP);
  Yimpm = zeros(Nfp*Nfaces,K);
  Yimpm(:) = 1 ./ Zimpm(:);
  Yimpp = zeros(Nfp*Nfaces,K);
  Yimpp(:) = 1 ./ Zimpp(:);
%
%  Homogeneous boundary conditions, Ez=0.
%
  Ebc = -E(vmapB);
  dE(mapB) = E(vmapB) - Ebc;
  Hbc = H(vmapB);
  dH(mapB) = H(vmapB) - Hbc;
%
%  Evaluate upwind fluxes.
%
  fluxE = 1 ./ (Zimpm + Zimpp) .* (nx.*Zimpp.*dH - dE);
  fluxH = 1 ./ (Yimpm + Yimpp) .* (nx.*Yimpp.*dE - dH);
%
%  Compute right hand sides of the PDE’s.
%
  rhsE = (-rx.*(Dr*H) + LIFT*(Fscale.*fluxE)) ./ eps;
  rhsH = (-rx.*(Dr*E) + LIFT*(Fscale.*fluxH)) ./ mu;

  return
end
