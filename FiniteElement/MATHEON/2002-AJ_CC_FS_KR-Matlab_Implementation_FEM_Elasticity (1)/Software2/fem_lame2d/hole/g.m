function sforce = g(x,n)
%G   Data on the Neumann boundary
%   Y = G(X,NORMAL) returns values of the normal-derivative at N discrete 
%   points on the Neumann boundary. This input data has to be choosen
%   by the user. X has dimension N x 2, NORMAL has dimension N x 2
%   containing the normal direction of the boundary at the corresponding
%   point, and Y has dimension N x 2.
%
%
%   See also FEM_LAME2D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <g.m> in $(HOME)/acfk/fem_lame2d/hole/

sforce = zeros(size(x,1),2);
sforce(find(n(:,2)==1),2) = 1;
