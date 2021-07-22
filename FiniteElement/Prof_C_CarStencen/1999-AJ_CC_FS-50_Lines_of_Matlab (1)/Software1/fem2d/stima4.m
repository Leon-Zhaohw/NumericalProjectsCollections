function M = stima4(vertices)
%STIMA4   Computes element stiffness matrix for parallelograms.
%   M = STIMA4(X) computes element stiffness matrix for parallelograms
%   analytically. X has dimension 4 x 2, where the first column gives 
%   the x-coordinates of the vertices and the second column their 
%   y-coordinates. M has dimension 4 x 4. The vertices are numbered 
%   anti-clockwise. 
%
%   This routine should not be modified.
%
%
%   See also FEM2D, STIMA3.

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <stima4.m> in $(HOME)/acf/fem2d/

D_Phi = [vertices(2,:)-vertices(1,:); vertices(4,:)-vertices(1,:)]';
B = inv(D_Phi'*D_Phi);
C1 = [2,-2;-2,2]*B(1,1)+[3,0;0,-3]*B(1,2)+[2,1;1,2]*B(2,2);
C2 = [-1,1;1,-1]*B(1,1)+[-3,0;0,3]*B(1,2)+[-1,-2;-2,-1]*B(2,2);
M = det(D_Phi) * [C1 C2; C2 C1] / 6;
