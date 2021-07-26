function [ Nv, VX, K, EToV ] = MeshGen1D ( xmin, xmax, K )

%*****************************************************************************80
%
%% MESHGEN1D generates a simple equidistant grid with K elements.
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
%  NV is the number of vertices.
%
  Nv = K + 1; 
%
%  VX contains the equally spaced vertex coordinates.
%
  VX = (1:Nv);
  for i = 1 : Nv
    VX(i) = ( xmax - xmin ) * ( i - 1 ) / ( Nv - 1 ) + xmin;
  end
%
%  EtoV is the Element-to-Vertex connectivity array.
%
  EToV = zeros ( K, 2 );

  for k = 1 : K
    EToV(k,1) = k; 
    EToV(k,2) = k + 1;
  end

  return
end
