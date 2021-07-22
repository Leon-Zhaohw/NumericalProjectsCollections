function [W,M] = u_d(x,lambda,mu)
%U_D   Data on the Dirichlet boundary
%   [W,M] = U_D(X,LAMBDA,MU) returns at N discrete points on the Dirichlet
%   boundary the direction for which the dispacement is given and the
%   corresponding values. If U is the displacement vector the Dirichlet
%   boundary condition is given by M*U = W. This input data has to be
%   choosen by the user. X has dimension N x 2, W has dimension 2*N x 1, and
%   M has dimension 2*N x 2. For this particular example, the Lame constants
%   LAMBDA and MU are also needed for the Dirichlet boundary condition.
%
%
%   See also FEM_LAME2D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <u_d.m> in $(HOME)/acfk/fem_lame2d/lshape_p1/ and
%                    $(HOME)/acfk/fem_lame2d/lshape_q1/

M = zeros(2*size(x,1),2);
W = zeros(2*size(x,1),1);
M(1:2:end,1) = 1;
M(2:2:end,2) = 1;
value = u_value(x,lambda,mu);
W(1:2:end,1) = value(:,1);
W(2:2:end,1) = value(:,2);
