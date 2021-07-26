% bump_demo.m
% This is a sample program illustrating how to:
%   - Generate a 3D planar triangle mesh, regular or irregular grid
%   - Determine interior and rows into exterior 
%   - Build and factor matrices for biharmonic and triharmonic systems with various
%     options
%   - Solve systems for new positions of mesh, given boundary conditions of
%     exterior
%
%   Alec Jacobson, Oct 29, 2010 
%   email with comments or questions: jacobson%cs.nyu.edu
%
%
% See corresponding paper: "Mixed finite elements for variational surface
% modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
% 2010
%
% Copyright 2010, Alec Jacobson, NYU
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for this sample program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Controls regularity of positions of points in domain:
%   'regular'   perfectly regular positions and connectivity
%   'irregular' parametrized variation from regular positions to completely
%               random
mesh_type = 'regular';
% Defines mass matrix type:
%   'barycentric' uses barycentric area
%   'voronoi' uses hybrid of voronoi area (for acute triangles) and constant 
%             area (for obtuse triangles)
masstype = 'voronoi';
% Determines boundary type for biharmonic system
%   'ext' region conditions (need to fix two rows into exterior)
%   'deriv' curve conditions (need to fix one row and specify tangents)
bi_bndtype = 'ext';
% Determines boundary type for triharmonic system
%   'extx' region conditions (need to fix three rows into exterior)
%   'extxy' extension curve conditions (need to fix two rows and specify bezier
%           curves at second row), experimental
%   'extxy-single' extension curve conditions (need to fix one row and specify 
%                  bezier curves), not implemented here
tri_bndtype = 'extx';
% Reduce system to only solve for u (eliminate auxillary variables by inverting
% the diagonal mass matrix) 
%   'flatten' reduce system, only valid for region conditions 
%             [Botsch and Kobbelt, 2004]
%   'no_flatten' do NOT reduce system, valid for all boundary conditions
reduction = 'no_flatten';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a 3D planar triangle mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(mesh_type,'regular'))
  % Regular triangle mesh grid with xRes, yRes points on each side
  % total number of vertices will be xRes * yRes
  % total number of triangle faces will be (xRes-1)*(yRes-1)*2
  xRes = 32;
  yRes = 32;
  % here we tell the function that we want a disc topology
  wrap = 0;
  % Output:
  %  F: list of triangles, each row is 3 indices to vertex list, thus #triangles x 3 array
  %  V: list of vertex positions #vertices x 2 array, (2D so far)
  %  res: legacy value
  %  edge_norms: list of norms for each edge, ordered like F for edges opposite
  %              vertex index, thus #triangles x 3 array
  [F,V,res,edge_norms] = create_regular_grid(xRes,yRes,wrap,wrap);

elseif(strcmp(mesh_type,'irregular'))
  % Irregular triangle mesh grid with xRes, yRes points on each side
  % These resolution parameters do not strictly govern how many points end up
  % in the resulting mesh, as the constrained delaunay triangulation may need
  % to add many steiner points. These are a rough control of the mesh
  % resolution.
  xRes = 32;
  yRes = 32;
  % here we tell the function that we want a disc topology
  wrap = 0;
  % How many "darts" to throw at random for each cell in regular grid of size
  % xRes by yRes
  points_per_cell = 1;
  % Minimum angle in degrees for angles on generated triangles:
  %   higher -> more extra points, better mesh quality
  %   lower -> fewer extra points, worse mesh quality
  min_angle = 30;
  % Variation from regular, how far from the regular positions we are allowed
  % to throw the "darts":
  %   higher -> more random, more irregular mesh
  %   lower -> closer to regular positions, regular mesh 
  dart_threshold = 0.5;

  % ( You must install Triangle and tell execute_triangle.m where your Triangle
  % program is)
  % Output:
  %  F: list of triangles, each row is 3 indices to vertex list, thus #triangles x 3 array
  %  V: list of vertex positions #vertices x 2 array, (2D so far)
  %  res: legacy value
  %  edge_norms: list of norms for each edge, ordered like F for edges opposite
  %              vertex index, thus #triangles x 3 array
  [F,V,res,edge_norms] = create_irregular_grid_with_min_angle( ...
    xRes, ... 
    yRes, ...
    points_per_cell, ...
    wrap, ...
    wrap, ...
    min_angle, ...
    false, ... % legacy parameter, keep as false
    dart_threshold);
end

% Define z-coordinates to be 0, now V is 3D
V = [V(:,1), V(:,2), zeros(size(V,1),1)];

% Display the input domain we've just generated
subplot_handle = subplot(2,2,1);
% clear subplot just in case
cla(subplot_handle);
tsurf(F,V);
title('Input domain');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine exterior and interior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the interior ("free" region) to be an annulus centered at (0.5,0.5) 
% R: major radius, r: minor radius
R = 0.4;
r = 0.2;
% Replace this with whatever should be the free region

% indices of all the vertices
indices = 1:size(V,1);
interior = indices( ...
  (V(:,1)-0.5).^2 + (V(:,2)-0.5).^2 < R^2 & ...
  (V(:,1)-0.5).^2 + (V(:,2)-0.5).^2 > r^2 );
exterior = indices( ...
  (V(:,1)-0.5).^2 + (V(:,2)-0.5).^2 >= R^2 | ...
  (V(:,1)-0.5).^2 + (V(:,2)-0.5).^2 <= r^2 );


% Biharmonic region system need 2 rows into the exterior, triharmonic 3 rows
% Define,
% Omega: interior from now on
% N0: First row into exterior, i.e. boundary between interior and exterior
% N1: 2nd row into exterior, i.e. boundary between Omega ∪ N0 and exterior \ N0
% N2: 3nd row into exterior, i.e. boundary between Omega ∪ N0 ∪ N1 and
%     exterior \ N0 \ N1
% (handle is a legacy term for exterior)
[Omega, N0, N1, N2, outside_region_of_interest] = ...
  layers_from_handle(size(V,1), F, exterior);

% Display exterior and interior on domain
subplot_handle = subplot(2,2,2);
% clear cubplot just in case
cla(subplot_handle);
display_domain(F,V,Omega,N0,N1,N2,outside_region_of_interest);
title('Paritioned domain');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build and factor matrices for bi-/tri-harmonic systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build and factor the biharmonic system matrix
% bi_L,bi_U,bi_P,bi_Q,bi_R are the factored parts of the system
% bi_S is the cotangent matrix, to be used in computing the rhs
% bi_M is the mass matrix, to be used in computing the rhs
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system( ...
  V, ...
  F, ...
  bi_bndtype, ...
  masstype, ...
  reduction, ...
  Omega, ...
  N0, ...
  N1);

% Build and factor the triharmonic system matrix
% tri_L,tri_U,tri_P,tri_Q,tri_R are the factored parts of the system
% tri_S is the cotangent matrix, to be used in computing the rhs
% tri_M is the mass matrix, to be used in computing the rhs
[tri_L,tri_U,tri_P,tri_Q,tri_R,tri_S,tri_M] = triharm_factor_system( ...
  V, ...
  F, ...
  tri_bndtype, ...
  masstype, ...
  reduction, ...
  Omega, ...
  N0, ...
  N1, ...
  N2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save original rest positions
rest_V = V;

% First move part of the exterior so that the solution is interesting

% Define the deformation to be raising the center disc up, depending on
% boundary conditions use tangent/bezier controls to specify noticable
% derivatives at the inner boundary

% Replace this with whatever should happen to the free region
V(exterior,:) = [ V(exterior,1) ...
  V(exterior,2) ...
  0.25*( (V(exterior,1)-0.5).^2 + (V(exterior,2)-0.5).^2 <= r^2)];

% Mimic a user interface where bezier curves have been specified at boundary
% vertices. Where BZ1 and BZ2 are #vertices by 3 arrays of 3D vectors defining
% a bezier curve in 3D at the corresponding vertex in V
%
%      BZ1
%   p0----->p1
%            |
%            | BZ2
%            ↓
%            p2
%                    
BZ1 = zeros(size(V,1),3);
BZ2 = zeros(size(V,1),3);

% Define tangent input
if(strcmp(bi_bndtype,'deriv'))
  % grab just inner ring boundary
  inner_N0 = N0( (V(N0,1)-0.5).^2 + (V(N0,2)-0.5).^2 <= r^2);
  % grad just outer ring boundary
  outer_N0 = N0( (V(N0,1)-0.5).^2 + (V(N0,2)-0.5).^2 >= R^2);

  % magnitude of tangents, roughly propotional to edge length
  inner_magnitude = 2*(1/xRes*2);
  outer_magnitude = (1/xRes*2);

  % Inner ring
  % point tangents away from center
  BZ1(inner_N0,:) = ...
    [(V(inner_N0,1)-0.5) (V(inner_N0,2)-0.5) 0.0*V(inner_N0,3)];
  % normalize and multiply by magnitude
  BZ1(inner_N0,:) = ...
    BZ1(inner_N0,:)./repmat(sqrt(sum(BZ1(inner_N0,:).^2,2)),1,3)*inner_magnitude;

  % Outer ring
  % point tangents toward from center
  BZ1(outer_N0,:) = ...
    [(0.5-V(outer_N0,1)) (0.5-V(outer_N0,2)) 0.0*V(outer_N0,3)];
  % normalize and multiply by magnitude
  BZ1(outer_N0,:) = ...
    BZ1(outer_N0,:)./repmat(sqrt(sum(BZ1(outer_N0,:).^2,2)),1,3)*outer_magnitude;
end

% Define bezier input
if(strcmp(tri_bndtype,'extxy'))
  % grab just inner ring boundary
  inner_N1 = N1( (V(N1,1)-0.5).^2 + (V(N1,2)-0.5).^2 <= r^2);
  % grad just outer ring boundary
  outer_N1 = N1( (V(N1,1)-0.5).^2 + (V(N1,2)-0.5).^2 >= R^2);

  % magnitude of the bezier segments, roughly propotional to edge length
  inner_magnitude = (1/xRes*2);
  outer_magnitude = (1/xRes*2)^2;

  % Inner ring
  % point first segment away from center
  BZ1(inner_N1,:) = ...
    [(V(inner_N1,1)-0.5) (V(inner_N1,2)-0.5) 0.0*V(inner_N1,3)];
  % normalize and multiply by magnitude
  BZ1(inner_N1,:) = ...
    BZ1(inner_N1,:)./repmat(sqrt(sum(BZ1(inner_N1,:).^2,2)),1,3)*inner_magnitude;
  % point second segement away from center and up
  BZ2(inner_N1,:) = ...
    [(V(inner_N1,1)-0.5) (V(inner_N1,2)-0.5) 1.0*V(inner_N1,3)];
  % normalize and multiply by magnitude
  BZ2(inner_N1,:) = ...
    BZ2(inner_N1,:)./repmat(sqrt(sum(BZ2(inner_N1,:).^2,2)),1,3)*inner_magnitude;

  % Outer ring
  % point first segment toward center
  BZ1(outer_N1,:) = ...
    [(0.5-V(outer_N1,1)) (0.5-V(outer_N1,2)) 0.0*V(outer_N1,3)];
  % normalize and multiply by magnitude
  BZ1(outer_N1,:) = ...
    BZ1(outer_N1,:)./repmat(sqrt(sum(BZ1(outer_N1,:).^2,2)),1,3)*magnitude;
  BZ2(outer_N1,:) = ...
    [(0.5-V(outer_N1,1)) (0.5-V(outer_N1,2)) 0.0*V(outer_N1,3)];
  % normalize and multiply by magnitude
  BZ2(outer_N1,:) = ...
    BZ2(outer_N1,:)./repmat(sqrt(sum(BZ2(outer_N1,:).^2,2)),1,3)*magnitude;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve systems for new positions of interior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate a space for the solutions to biharmonic and triharmonic systems
bi_V = zeros(size(V));
tri_V = zeros(size(V));
% fill in exterior with known "fixed" values
bi_V(exterior,:) = V(exterior,:);
tri_V(exterior,:) = V(exterior,:);

% Build RHS and use factor system matrix to solve for new positions in each
% coordinate
bi_V = biharm_solve_with_factor( ...
  bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
  F, V, Omega, N0, N1, bi_bndtype, reduction,BZ1,rest_V);

% Display solution
subplot_handle = subplot(2,2,3);
% clear cubplot just in case
cla(subplot_handle);
display_domain(F,bi_V,Omega,N0,N1,N2,outside_region_of_interest,1);
view(3)
zoom out
zoom(1.5)
title('Biharmonic solution');
axis equal; 
axis vis3d;

% Build RHS and use factor system matrix to solve for new positions in each
% coordinate
tri_V = triharm_solve_with_factor( ...
  tri_L, tri_U, tri_P, tri_Q, tri_R, tri_S, tri_M, ...
  F, V, Omega, N0, N1, N2, tri_bndtype, reduction,BZ1, BZ2, rest_V);

% Display solution
subplot_handle = subplot(2,2,4);
% clear cubplot just in case
cla(subplot_handle);
display_domain(F,tri_V,Omega,N0,N1,N2,outside_region_of_interest,1);
view(3)
zoom out
zoom(1.5)
title('Triharmonic solution');
axis equal; 
axis vis3d;
