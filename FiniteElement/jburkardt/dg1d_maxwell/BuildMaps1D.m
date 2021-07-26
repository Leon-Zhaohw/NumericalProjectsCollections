function [ vmapM, vmapP, vmapB, mapB ] = BuildMaps1D ( )

%*****************************************************************************80
%
%% BUILDMAPS1D creates connectivity and boundary tables for nodes.
%
%  Discussion:
%
%    For each of the K elements there are nodes associated with N+1 degrees
%    of freedom.
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
%  number volume nodes consecutively
%
  nodeids = reshape(1:K*Np, Np, K);
  vmapM   = zeros(Nfp, Nfaces, K); 
  vmapP   = zeros(Nfp, Nfaces, K); 

  for k1=1:K
    for f1=1:Nfaces
%
%  find index of face nodes with respect to volume node ordering
% 
      vmapM(:,f1,k1) = nodeids(Fmask(:,f1), k1);
    end
  end

  for k1=1:K
    for f1=1:Nfaces
%
%  Find neighbor
%
      k2 = EToE(k1,f1);
      f2 = EToF(k1,f1);
%   
%  Find volume node numbers of left and right nodes 
%
      vidM = vmapM(:,f1,k1);
      vidP = vmapM(:,f2,k2);
    
      x1  = x(vidM);
      x2  = x(vidP);
%   
%  Compute distance matrix
%
      D = (x1 -x2 ).^2;
      if (D<NODETOL)
        vmapP(:,f1,k1) = vidP;
      end;
    end
  end

  vmapP = vmapP(:);
  vmapM = vmapM(:);
%
%  Create list of boundary nodes
%
  mapB = find(vmapP==vmapM);
  vmapB = vmapM(mapB);
%
%  Create specific left (inflow) and right (outflow) maps
%
  mapI = 1;
  mapO = K*Nfaces;
  vmapI = 1;
  vmapO = K*Np;

  return
end

