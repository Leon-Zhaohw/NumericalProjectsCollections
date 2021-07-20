function volforce = f(x);
%F   Volume force in considered domain.
%   Y = F(X) returns values of forces at N discrete points in the considered
%   domain. This input data has to be chosen by the user. X has dimension N
%   x 2 and Y has dimension N x 2.
%
%
%   See also FEM_LAME2D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <f.m> in $(HOME)/acfk/fem_lame2d/cooks/

volforce = zeros(size(x,1),2);
