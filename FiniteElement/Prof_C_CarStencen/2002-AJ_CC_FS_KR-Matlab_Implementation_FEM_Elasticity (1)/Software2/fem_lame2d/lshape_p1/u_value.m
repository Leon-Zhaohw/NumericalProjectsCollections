function value = u_value(x,lambda,mu) 
%U_VALUE   Known solution of the L-Shape problem.
%   U_VALUE(X,LAMBDA,MU) returns the known value of the L-shape problem
%   at the positions given in X, a N x 2 matrix. LAMBDA and MU are the
%   Lame constants.
%
%
%   See also FEM_LAME2D and U_D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <u_value.m> in $(HOME)/acfk/fem_lame2d/lshape_p1/ and
%                        $(HOME)/acfk/fem_lame2d/lshape_q1/

[phi,r] = cart2pol(x(:,1),x(:,2));
alph = .544483737;
omega = 3*pi/4;
C_1 = -cos((alph+1)*omega)/cos((alph-1)*omega);
C_2 = 2*(lambda+2*mu)/(lambda+mu);
ut = (1/(2*mu)) * r.^alph .*((alph+1)*sin((alph+1)*phi)+ ...
    (C_2+alph-1)*C_1*sin((alph-1)*phi));
ur = (1/(2*mu))*r.^alph .* (-(alph+1)*cos((alph+1)*phi)+ ...
    (C_2-(alph+1))*C_1*cos((alph-1)*phi));
value = [ur .* cos(phi) - ut .* sin(phi), ur .* sin(phi) + ut .* cos(phi)];
