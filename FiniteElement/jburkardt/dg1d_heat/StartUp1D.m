%*****************************************************************************80
%
%% STARTUP1D is a script that sets building operators, grid, metric and connectivity. 
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
%  Definition of global constants
%
Globals1D; 

NODETOL = 1e-10;
Np = N + 1; 
Nfp = 1; 
Nfaces = 2;
%
%  Compute basic Legendre-Gauss-Lobatto grid
%
r = JacobiGL ( 0, 0, N );
%
%  Build reference element matrices
%
V  = Vandermonde1D ( N, r ); 
invV = inv ( V );
Dr = Dmatrix1D ( N, r, V );
%
%  Create surface integral terms
%
LIFT = Lift1D ( );
%
%  Build coordinates of all the nodes
%
va = EToV(:,1)'; 
vb = EToV(:,2)';
x = ones(N+1,1) * VX(va) + 0.5 * (r+1) * (VX(vb)-VX(va));
%
%  Calculate geometric factors
%
[rx,J] = GeometricFactors1D ( x, Dr );
%
%  Compute masks for edge nodes
%
fmask1 = find( abs(r+1) < NODETOL)'; 
fmask2 = find( abs(r-1) < NODETOL)';
Fmask  = [fmask1;fmask2]';
Fx = x(Fmask(:), :);
%
%  Build surface normals and inverse metric at surface
%
[nx] = Normals1D();
Fscale = 1 ./ (J(Fmask,:));
%
%  Build connectivity matrix
%
[EToE, EToF] = Connect1D ( EToV );
%
%  Build connectivity maps
%
[ vmapM, vmapP, vmapB, mapB ] = BuildMaps1D ( );

