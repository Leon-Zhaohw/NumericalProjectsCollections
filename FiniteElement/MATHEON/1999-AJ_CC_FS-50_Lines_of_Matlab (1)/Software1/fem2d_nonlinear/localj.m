function b = localj(vertices,U)
% LOCALJ   Computes local functional J(U,V) in considered problem.
%    B = LOCALJ(VERTICES,U) computes the discretized local functional J(U,V_j)
%    over the triangle specified by VERTICES for given U and for V_j being the
%    three local hat functions at the vertices of the triangle. This M-file
%    is problem dependent and must be chosen by the user. VERTICES has
%    dimension 3 X 2 where the first column gives the x-coordinates of the 3
%    vertices and the second column their y-coordinates. U is a 3 x 1
%    dimensional array consisting of the values at the corresponding
%    vertices. B is a 3 x 1 dimensional array returning the value of the
%    functional J(U,V_j) for the three corresponding vertices.
%
%
%    See also FEM2D_NONLINEAR and LOCALDJ.

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <localj.m> in $(HOME)/acf/fem2d_nonlinear/

Eps = 1/100;
G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
Area = det([ones(1,3);vertices']) / 2;
b=Area*((Eps*G*G'-[2,1,1;1,2,1;1,1,2]/12)*U+ ...
    [4*U(1)^3+ U(2)^3+U(3)^3+3*U(1)^2*(U(2)+U(3))+2*U(1) ...
      *(U(2)^2+U(3)^2)+U(2)*U(3)*(U(2)+U(3))+2*U(1)*U(2)*U(3);
  4*U(2)^3+ U(1)^3+U(3)^3+3*U(2)^2*(U(1)+U(3))+2*U(2) ...
      *(U(1)^2+U(3)^2)+U(1)*U(3)*(U(1)+U(3))+2*U(1)*U(2)*U(3);
  4*U(3)^3+ U(2)^3+U(1)^3+3*U(3)^2*(U(2)+U(1))+2*U(3) ...
      *(U(2)^2+U(1)^2)+U(2)*U(1)*(U(2)+U(1))+2*U(1)*U(2)*U(3)]/60);
