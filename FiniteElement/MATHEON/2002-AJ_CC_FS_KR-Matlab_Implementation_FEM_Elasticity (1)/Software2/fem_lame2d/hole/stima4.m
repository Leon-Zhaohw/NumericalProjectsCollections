function stima4 = stima4(vertices,lambda,mu)
%STIMA4   Computes element stiffness matrix for parallelograms.
%   M = STIMA4(X,LAMBDA,MU) computes element stiffness matrix for
%   parallelograms. X has dimension 4 x 2, where the first column gives the
%   x-coordinates of the vertices and the second column their
%   y-coordinates.The vertices are numbered anti-clockwise. LAMBDA and MU
%   are the Lame constants. M has dimension 8 x 8.
%
%   This routine should not be modified.
%
%
%   See also FEM_LAME2D and STIMA3.

%    J. Alberty, C. Carstensen and S. A. Funken  07-03-00
%    File <stima4.m> in $(HOME)/acfk/fem_lame2d/cooks/ and
%                       $(HOME)/acfk/fem_lame2d/lshape_p1/ and
%                       $(HOME)/acfk/fem_lame2d/lshape_q1/ and
%                       $(HOME)/acfk/fem_lame2d/hole/

R_11 = [2 -2 -1  1; -2  2  1 -1; -1  1 2 -2;  1 -1 -2  2]/6;
R_12 = [1  1 -1 -1; -1 -1  1  1; -1 -1 1  1;  1  1 -1 -1]/4;
R_22 = [2  1 -1 -2;  1  2 -2 -1; -1 -2 2  1; -2 -1  1  2]/6;
F = inv([vertices(2,:)-vertices(1,:); vertices(4,:)-vertices(1,:)]);
L = [lambda+2*mu,lambda,mu];
stima4 = zeros(8,8); 
E = F'*[L(1),0;0,L(3)]*F;
stima4(1:2:8,1:2:8) = E(1,1)*R_11 +E(1,2)*R_12 +E(2,1)*R_12' +E(2,2)*R_22; 
E = F'*[L(3),0;0,L(1)]*F;
stima4(2:2:8,2:2:8) = E(1,1)*R_11 +E(1,2)*R_12 +E(2,1)*R_12' +E(2,2)*R_22; 
E = F'*[0,L(3);L(2),0]*F;
stima4(2:2:8,1:2:8) = E(1,1)*R_11 +E(1,2)*R_12 +E(2,1)*R_12' +E(2,2)*R_22; 
stima4(1:2:8,2:2:8) = stima4(2:2:8,1:2:8)';
stima4 = stima4/det(F);
