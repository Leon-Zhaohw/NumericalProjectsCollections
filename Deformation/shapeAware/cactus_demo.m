% This is a script that demos computing automatic, monotonic weights for point
% handles on a 2D shape, according to "Smooth Shape-Aware Functions with
% Controlled Extrema" by Alec Jacobson, Tino Weinkauf, and Olga Sorkine 2012
%
% This file and any included files (unless otherwise noted) are copyright Alec
% Jacobson. Email jacobson@inf.ethz.ch if you have questions
%
% Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
%

% NOTE: Please contact Alec Jacobson, jacobson@inf.ethz.ch before
% using this code outside of an informal setting, i.e. for comparisons.

% Close all open figure windows
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load a mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input mesh source: *.obj, *.off, *.poly, or *.png
mesh_source = 'cactus.obj';
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
fprintf( ...
  ['\nCLICK on mesh at each location where you would like to add a ' ...
  'point handle.\n' ...
  'Press ENTER when finished.\n\n']);
% User clicks many times on mesh at locations of control points
C = get_control_points(V,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bind controls to mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: This computes the "normalized" or "optimized" version of our monotonic
% weights, *not* the full solution which solve for all weights simultaneously
% and enforce partition of unity as a proper contstraint. 

% Compute boundary conditions
[b,bc] = boundary_conditions(V,F,C);
% Compute [Botsch & Kobbelt 2004] weights
Wb = harmonic(V,F,b,bc,2);
% Compute [Jacobson et al. 2011] weights
Wbbw = biharmonic_bounded(V,F,b,bc);
% Normalize weights
Wbbw = bsxfun(@rdivide,Wbbw,sum(Wbbw,2));
% Compute weights
Wm = monotonic_biharmonic(V,F,b,bc);
% Normalize weights
Wm = bsxfun(@rdivide,Wm,sum(Wm,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deform mesh via controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of isointervals to display
nin = 20;

% Display mesh and control points and allow user to interactively deform mesh
% and view weight visualizations

% NOTE: Deformation uses Dual Quaternion Skinning rather than Linear Blend
% Skinning to avoid confusing scaling artifacts due to Linear Blending
% rotations regardless of weights.

figure(1);
simple_deform(V,F,C,Wb,'InterpMode','DQLBS','ShowWeightVisualization');
set(gca, 'visible', 'off');
title('Unconstrained Biharmonic','Visible','on');
% Compute a pseudo color "heat" map to visualize weights
[HM,HC] = weights_colormap(Wb(:),nin);
colormap(HM);
caxis(HC);
colorbar;

figure(2);
simple_deform(V,F,C,Wbbw,'InterpMode','DQLBS','ShowWeightVisualization');
set(gca, 'visible', 'off');
title('Bounded Biharmonic','Visible','on');
[HM,HC] = weights_colormap(Wbbw(:),nin);
colormap(HM);
caxis(HC);
colorbar;

figure(3);
simple_deform(V,F,C,Wm,'InterpMode','DQLBS','ShowWeightVisualization');
set(gca, 'visible', 'off');
title('Monotonic Biharmonic','Visible','on');
[HM,HC] = weights_colormap(Wm(:),nin);
colormap(HM);
caxis(HC);
colorbar;

% Tile figures across screen
tilefigs([1 inf]);
