function stima3=stima3(vertices,lambda,mu)
%STIMA3   Computes element stiffness matrix for triangles.
%   M = STIMA3(X,LAMBDA,MU) computes element stiffness matrix for
%   triangles. The coordinates of the vertices are stored in X. LAMBDA
%   and MU are the Lame constants.
%
%   This routine should not be modified.
%
%
%   See also FEM_LAME2D and STIMA4.

%    J. Alberty, C. Carstensen and S. A. Funken  07-03-00
%    File <stima3.m> in $(HOME)/acfk/fem_lame2d/cooks/ and
%                       $(HOME)/acfk/fem_lame2d/lshape_p1/ and
%                       $(HOME)/acfk/fem_lame2d/lshape_q1/ and
%                       $(HOME)/acfk/fem_lame2d/hole/

PhiGrad = [1,1,1;vertices']\[zeros(1,2);eye(2)];
R = zeros(3,6);
R([1,3],[1,3,5]) = PhiGrad';
R([3,2],[2,4,6]) = PhiGrad';
C = mu*[2,0,0;0,2,0;0,0,1] +lambda*[1,1,0;1,1,0;0,0,0];
stima3 = det([1,1,1;vertices'])/2*R'*C*R;
