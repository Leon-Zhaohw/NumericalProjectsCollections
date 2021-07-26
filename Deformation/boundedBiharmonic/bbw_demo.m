% This is a script that demos computing Bounded Biharmonic Weights
% automatically for a 2D shape.
%
% This file and any included files (unless otherwise noted) are copyright Alec
% Jacobson. Email jacobson@inf.ethz.ch if you have questions
%
% Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
%

% NOTE: Please contact Alec Jacobson, jacobson@inf.ethz.ch before
% using this code outside of an informal setting, i.e. for comparisons.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load a mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input mesh source: *.obj, *.off, *.poly, or *.png
mesh_source = 'woody.obj';
% should input mesh be upsampled
upsample_mesh = false;

if(~isempty(regexp(mesh_source,'\.(off|obj)$')))
  % load a mesh from an OBJ
  [V,F] = load_mesh(mesh_source);
  % only keep x and y coordinates, since we're working only in 2D
  V = V(:,1:2);
elseif ~isempty(regexp(mesh_source,'\.poly$'))
  % load a mesh from a .POLY polygon file format
  % Triangulate in two-passes. First pass with just angle constraint forces
  % triangles near the boundary to be small, but internal triangles will be very
  % graded
  [V,F] = triangle(mesh_source,'Quality',30);
  % phony z-coordinate
  V = [V, zeros(size(V,1),1)];
  % compute minimum angle 
  min_area = min(doublearea(V,F))/2;
  % Use minimum area of first pass as maximum area constraint of second pass for
  % a more uniform triangulation. probably there exists a good heuristic for a
  % maximum area based on the input edge lengths, but for now this is easy
  % enough
  [V,F] = triangle(mesh_source,'Quality',30,'MaxArea',min_area);
elseif ~isempty(regexp(mesh_source,'\.png$'))
  % load a mesh from a PNG image with transparency
  [V,F] = png2mesh(mesh_source,1,50);
end

% upsample each triangle
if(upsample_mesh)
  [V,F] = upsample(V,F);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Place controls on mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% display mesh
tsurf(F,V)
axis equal;
fprintf( ...
  ['\nCLICK on mesh at each location where you would like to add a ' ...
  'point handle.\n' ...
  'Press ENTER when finished.\n\n']);
% User clicks many times on mesh at locations of control points
try
  [Cx,Cy] = getpts;
catch e
  % quit early, stop script
  return
end
% store control points in single #P by 2 list of points
C = [Cx,Cy];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bind controls to mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: This computes the "normalized" or "optimized" version of BBW, *not* the
% full solution which solve for all weights simultaneously and enforce
% partition of unity as a proper contstraint. 

% Compute boundary conditions
[b,bc] = boundary_conditions(V,F,C);
% Compute weights
if(exist('mosekopt','file'))
  % if mosek is installed this is the fastest option
  W = biharmonic_bounded(V,F,b,bc,'conic');
else
  % else this uses the default matlab quadratic programming solver
  W = biharmonic_bounded(V,F,b,bc,'quad');
end
% Normalize weights
W = W./repmat(sum(W,2),1,size(W,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deform mesh via controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display mesh and control points and allow user to interactively deform mesh
% and view weight visualizations

% Commented out are examples of how to call "simple_deform.m" with various options
simple_deform(V,F,C,W,'ShowWeightVisualization');
% interactively deform point controls
%simple_deform(V,F,C,W)

