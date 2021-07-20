function [W,M] = u_d(x)
%U_D   Data on the Dirichlet boundary
%   [W,M] = U_D(X) returns at N discrete points on the Dirichlet boundary
%   the direction for which the dispacement is given and the corresponding
%   values. If U is the displacement vector the Dirichlet boundary condition
%   is given by M*U = W. This input data has to be choosen by the user. X
%   has dimension N x 2, W has dimension 2*N x 1, and M has dimension 2*N
%   x 2.
%
%
%   See also FEM_LAME2D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <u_d.m> in $(HOME)/acfk/fem_lame2d/hole/

M = zeros(2*size(x,1),2);
W = zeros(2*size(x,1),1);
% symmetry conditions on the x-axis
temp = find(x(:,1)>0 & x(:,2)==0);
M(2*temp-1,2) = 1;
% symmetry conditions on the y-axis
temp = find(x(:,2)>0 & x(:,1)==0); 
M(2*temp-1,1) = 1;
