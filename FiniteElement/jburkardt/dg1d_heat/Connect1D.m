function [ EToE, EToF ] = Connect1D ( EToV )

%*****************************************************************************80
%
%% CONNECT1D builds global connectivity arrays for the 1D grid.
%
%  Discussion:
%
%    The calculation is based on the standard EToV input array from the
%    grid generator.
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
  Nfaces = 2;
%
%  Find number of elements and vertices
%
  K = size(EToV,1); 
  TotalFaces = Nfaces*K; 
  Nv = K+1;
%
%  List of local face to local vertex connections
%
  vn = [1,2];
%
%  Build global face to node sparse array
%
  SpFToV = spalloc(TotalFaces, Nv, 2*TotalFaces);
  sk = 1;
  for k=1:K
    for face=1:Nfaces
      SpFToV( sk, EToV(k, vn(face))) = 1;
      sk = sk+1;
    end
  end
%
%  Build global face to global face sparse array
%
  SpFToF = SpFToV*SpFToV' - speye(TotalFaces);
%
%  Find complete face to face connections
% 
  [faces1, faces2] = find(SpFToF==1);
%
%  Convert face global number to element and face numbers
%
  element1 = floor( (faces1-1)/Nfaces )  + 1;
  face1    =   mod( (faces1-1), Nfaces ) + 1;
  element2 = floor( (faces2-1)/Nfaces )  + 1;
  face2    =   mod( (faces2-1), Nfaces ) + 1;
%
%  Rearrange into Nelements x Nfaces sized arrays
%
  ind = sub2ind([K, Nfaces], element1, face1);
  EToE      = (1:K)'*ones(1,Nfaces);
  EToF      = ones(K,1)*(1:Nfaces);
  EToE(ind) = element2;
  EToF(ind) = face2;

  return
end

