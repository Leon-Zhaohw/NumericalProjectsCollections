function M = stima3(vertices)
%STIMA3   Computes element stiffness matrix for simplex.
%   M = STIMA3(X) computes element stiffness matrix for simplex, 
%   i.e. for triangles in two dimensions (d=2) and tetraeder in 
%   three dimensions (d=3). The coordinates of the vertices are stored 
%   in X. X has dimension (d+1) x d. In two-dimension, the vertices 
%   are numbered anti-clockwise. M has dimension (d+1) x (d+1). In 
%   three-dimension, the vertices are numbered s.t. max(eig(M)) > 0. 
%
%   This routine should not be modified.

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <stima3.m> in $(HOME)/acf/fem2d/ and
%                    in $(HOME)/acf/fem3d/ and
%                    in $(HOME)/acf/fem2d_heat/

d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
