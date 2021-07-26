function [ p, t, stats ] = mesh2d ( node, edge, hdata, options )

%*****************************************************************************80
%
%% mesh2d(): 2D unstructured mesh generation for a polygon.
%
%  Discussion:
%
%    A 2D unstructured triangular mesh is generated based on a piecewise-
%    linear geometry input.  The polygon can contain an arbitrary number of 
%    cavities.  An iterative method is implemented to optimise mesh quality. 
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 December 2020
%
%  Author:
%
%    Darren Engwirda
%
%  Input:
%
%    real NODE(N,2), geometry nodes, normally specified in consecutive order, 
%    such that NODE(2,:) is joined with NODE(1,:), and so on.
%
%    integer EDGE(NE,2), an optional input.  EDGE defines the connectivity 
%    between the points in NODE as a list of edges:
%      edge = [n1 n2; n2 n3; etc]
%    If EDGE is specified, it is not required that NODE be consecutive.
%    Otherwise, EDGE may be input as [].
%
%    HDATA
%
%    OPTIONS
%
%  Output:
%
%    real P(N,2), nodal XY co-ordinates.
%
%    real T(M,3), triangles as indices into P, defined with a
%    counter-clockwise node ordering.
%
%
% An element size function is automatically generated based on the 
% complexity of the geometry. Generally this produces meshes with the 
% fewest number of triangles.
%
% LONG SYNTAX:
%
%  [p,t] = mesh2d(node,edge,hdata,options);
%
% Blank arguments can be passed using the empty placeholder "[]".
%
% HDATA is a structure containing user defined element size information. 
% HDATA can include the following fields:
%
%  hdata.hmax  = h0;                   Max allowable global element size.
%  hdata.edgeh = [e1,h1; e2,h2; etc];  Element size on specified geometry 
%                                      edges.
%  hdata.fun   = 'fun' or @fun;        User defined size function.
%  hdata.args  = {arg1, arg2, etc};    Additional arguments for HDATA.FUN.
%
% Calls to user specified functions must accept vectorised input of the 
% form H = FUN(X,Y,ARGS{:}), where X,Y are the xy coordinates where the
% element size will be evaluated and ARGS are optional additional arguments 
% as passed by HDATA.ARGS.
%
% An automatic size function is always generated to ensure that the
% geometry is adequately resolved. The overall size function is the minimum
% of the user specified and automatic functions.
%
% OPTIONS is a structure array that allows some of the "tuning" parameters
% used in the solver to be modified:
%
%   options.mlim   : The convergence tolerance. The maximum percentage 
%                    change in edge length per iteration must be less than 
%                    MLIM { 0.02, 2.0% }. 
%   options.maxit  : The maximum allowable number of iterations { 20 }.
%   options.dhmax  : The maximum allowable (relative) gradient in the size 
%                    function { 0.3, 30.0% }.
%   options.output : Displays the mesh { TRUE }.
%   options.stats  : Displays the mesh statistics { FALSE }.
%
% STATS is an undocumented output used in debugging. Returns the algorithm 
% statistics usually printed to screen as a structure.
%
  if (nargin<4)
    options = [];
    options.output = true;
    options.stats = false;
    if (nargin<3)
      hdata = [];
      if (nargin<2)
        edge = [];
      end
    end
  end
%
%  Assume 1 face containing all edges
%
  [ p, t, junk, stats ] = meshfaces ( node, edge, [], hdata, options );

  return
end

