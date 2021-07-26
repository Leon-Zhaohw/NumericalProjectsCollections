% This is a script that demos computing Bounded Biharmonic Weights
% automatically for the 2D alligator as seen in figure 1 of "Bounded Biharmonic
% Weights for Real-time Deformation" by Jacobson et al.
%
% This file and any included files (unless otherwise noted) are copyright Alec
% Jacobson. Email jacobson@inf.ethz.ch if you have questions
%
% Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
%

% Zip up using:
% >> C = depends('alligator_demo');
% >> C = C(cellfun(@isempty,strfind(C,'opt/local/mosek')));
% >> C = ...
%   cat(1,C,{ ...
%     'README'; ...
%     'alligator.obj'; ...
%     'alligator.png'; ...
%     'alligator-skeleton-cage-points.tgf'});
% >> zip('alligator_demo.zip',C);
% >> fprintf('This package should contain\nalligator_demo/\n');
% >> N = regexprep(C,'^.*\/','');
% >> fprintf('  %s\n',N{:});

% NOTE: Please contact Alec Jacobson, jacobson@inf.ethz.ch before
% using this code outside of an informal setting, i.e. for comparisons.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load a alligator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input mesh source: *.obj, *.off, *.poly, or *.png
mesh_source = 'alligator.png';
image_source = 'alligator.png';

if(~isempty(regexp(mesh_source,'\.(off|obj)$')))
  % load a mesh from an OBJ
  [V,F] = load_mesh(mesh_source);
  % only keep x and y coordinates, since we're working only in 2D
  V = V(:,1:2);
elseif ~isempty(regexp(mesh_source,'\.png$'))
  % load a mesh from a PNG image with transparency
  if(exist('mosekopt','file'))
    [V,F] = png2mesh(mesh_source,1,500);
  else
    % matlab's solver is slow so use a coarser mesh
    [V,F] = png2mesh(mesh_source,80,300);
  end
end
tsurf(F,V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load controls for alligator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input handle source: *.tgf
handle_source = 'alligator-skeleton-cage-points.tgf';
% read in handles from .tgf file 
[C,E,P,BE,CE,PE] = readTGF(handle_source);
% only keep x and y
C = C(:,1:2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample point handles, bone edges and cages, then remesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
original_V = V;
original_F = F;
[V,F,OV,OE,OV1,OE1,epsilon] = remesh_at_handles(V,F,C,P,BE,CE);

% show a plot before cleaning up triangle input
subplot(3,1,1);
plot([OV1(OE1(:,1),1) OV1(OE1(:,2),1)]',[OV1(OE1(:,1),2) OV1(OE1(:,2),2)]','-','LineWidth',1);
hold on;
plot(OV1(:,1),OV1(:,2),'.');
hold off;
axis equal;
title('Outline + handle samples');
% show a plot of what's sent to triangle
subplot(3,1,2);
plot([OV(OE(:,1),1) OV(OE(:,2),1)]',[OV(OE(:,1),2) OV(OE(:,2),2)]','-','LineWidth',1);
hold on;
plot(OV(:,1),OV(:,2),'.');
hold off;
axis equal;
title('Input to triangle (collapsed points + split edges)');
% show a plot of the result
subplot(3,1,3);
tsurf(F,V);
hold on;
plot([OV(OE(:,1),1) OV(OE(:,2),1)]',[OV(OE(:,1),2) OV(OE(:,2),2)]','-','LineWidth',2);
hold off;
axis equal;
title('Output of triangle');
tilefigs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bind controls to mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: This computes the "normalized" or "optimized" version of BBW, *not* the
% full solution which solve for all weights simultaneously and enforce
% partition of unity as a proper contstraint. 

% Compute boundary conditions
[b,bc] = boundary_conditions(V,F,C,P,BE,CE);
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

% only keep faces of part of mesh whose centroids are in original shape
[B,L] = ordered_outline(original_F);
[xin,iin,cin] = faces_in_polygon(V,F,original_V(B,1),original_V(B,2));
% just limit faces so exterior faces don't show up
% would be better to remap vertices, faces and weights...

% Display mesh and control points and allow user to interactively deform mesh
% and view weight visualizations

% interactively deform point controls
figure;
global g_Deform;
gid = simple_deform(V,F(cin,:),C,W,'PointHandles',P,'BoneEdges',BE,'CageEdges',CE);
set(g_Deform(gid).tsh,'FaceColor',[0.1 0.5 0.1])
set(g_Deform(gid).tsh,'EdgeColor',[0.1 0.5 0.1])
tilefigs([inf,1]);

% EXPERIMENTAL 
%[im,map,alpha] = imread(image_source);
%% treat alpha as white
%im = im2double(im).*repmat(im2double(alpha),[1 1 size(im,3)]) + 1.*repmat(1-im2double(alpha),[1 1 size(im,3)]);
%% get height and width
%h = size(im,1);
%w = size(im,2);
%% size of grid cells
%step = 3;
%[X,Y] = meshgrid(1:step:w,1:step:h);
%% matlab likes to display warped images with backwards frame
%flipV = [V(:,1) h-V(:,2)];
%original_flip_V = [original_V(:,1) h-original_V(:,2)];
%% also flip control point positions
%flipC = [C(:,1) h-C(:,2)];
%% get interpolant for weight functions
%A = TriScatteredInterpVector(flipV,W);
%% interpolate weights at grid
%WXY = reshape(A(X,Y),prod(size(X)),size(W,2));
%
%% white-out NaNs in image
%[fullX,fullY] = meshgrid(1:w,1:h);
%nan_interp = TriScatteredInterp([X(:) Y(:)],double(any(isnan(WXY),2)));
%%im_nan = floor(repmat(nan_interp(fullX,fullY),[1 1 3]));
%%im_masked = im.*(1-im_nan) + 1*(im_nan);
%im_in = inpolygon(fullX, fullY, original_flip_V([B(L(1):(L(1+1)-1)) B(L(1))],1),original_flip_V([B(L(1):(L(1+1)-1)) B(L(1))],2));
%im_in = imerode(im_in,strel('disk',2));
%im_in = repmat(im_in,[1 1 3]);
%im_masked = im.*(im_in) + 1*(1-im_in);
%
%deform([],[],flipC,'DeformMethod','skinning',WXY,'PointHandles',P,'BoneEdges',BE,'CageEdges',CE,'Image',im_masked,X,Y);
