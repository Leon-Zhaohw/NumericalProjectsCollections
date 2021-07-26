%*****************************************************************************80
%
%% GLOBALS1D declares global variables
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
%  Dr(1:N+1,1:N+1) is the differentiation matrix, equaling inv(VR)*V,
%  where V is the Vandermonde matrix, and VR is the gradient of V.
%
global Dr
%
%  EToE(1:K,1:Nfaces) is the element-to-element connectivity array.
%
global EToE 
%
%  EtoF(1:K,1:Nfaces) is the element-to-face connectivity array.
%
global EToF
%
%  Fmask(left edges + right edges) is a mask for edge nodes, that is, nodes 
%  whose "r" value is -1 or +1.
%
global Fmask
%
%  Fscale(?) is the inverse metric evaluated at the surface.
%
global Fscale
%
%  Fx(left edges + right edges) is the list of edge nodes: Fx = x(Fmask(:),:).
%
global Fx
%
%  invV(1:N+1,1:N+1) is the inverse of V, the Vandermonde matrix.
%
global invV
%
%  J(?,?) is the metric element for the local mapping of the 1D elements.
%
global J
%
%  K is the number of elements.
%
global K
%
%  LIFT(1:N+1,1:Nfaces) is the surface integral term in the DG formulation.
%
global LIFT
%
%  mapB(?): list of boundary nodes.
%
global mapB
%
%  mapI(?): left (inflow) map.
%
global mapI
%
%  mapO(?): right (outflow) map.
%
global mapO
%
%  N is the order of the polynomials used for approximation.
%
global N
%
%  Nfaces is 2, the number of "faces" each element has.
%
global Nfaces
%
%  Nfp is set to 1.  WHAT DOES IT MEAN?
%
global Nfp
%
%  NODETOL is a very small tolerance for locating nodes.
%
global NODETOL
%
%  Np is N+1, the actual number of coefficients in a polynomial of degree N.
%
global Np
%
%  nx(1:2,K) is the surface normals.  
%  For this 1D problem, they are -1 and +1 for each element.
%
global nx
%
%  r(1:N+1) contains the Gauss-Lobatto quadrature points for [-1,+1] of order N.
%
global r
%
%  rk4a(1:5): Runge-Kutta coefficients:
%
global rk4a
%
%  rk4b(1:5): Runge-Kutta coefficients:
%
global rk4b
%
%  rk4c(1:5): Runge-Kutta coefficients:
%
global rk4c
%
%  rx(?,?) is the inverse of J.
%
global rx
%
%  V(1:N+1,1:N+1) is the Vandermonde matrix.
%
global V
global vmapB
global vmapI
%
%  vmapM(?): index of face nodes with respect to volume node ordering.
%
global vmapM
global vmapO
global vmapP
%
%  VX(1:K+1) contains the coordinates of the vertices, 
%  which are the endpoints of elements.
%
global VX
%
%  x(1:N+1,K) contains the coordinates of all the nodes,
%  indexed by local node index and element.
%
global x
%
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];

rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];

rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

