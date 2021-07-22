function DirichletRandwert = u_d(x)
%U_D   Data on the Dirichlet boundary
%   Y = U_D(X) returns function values at N discrete points  
%   on the Dirichlet boundary. This input data has to be choosen  
%   by the user. X has dimension N x 3 and Y has dimension N x 1.
%
%
%   See also FEM3D, F, G.
%

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <u_d.m> in $(HOME)/acf/fem3d/
%    This Dirichlet boundary data is used to compute Fig. 6 in 
%    "Remarks around "50 lines of Matlab: Short finite element 
%    implementation"

DirichletRandwert = -150*ones(size(x,1),1);
ind = find( x(:,1) >= 2.0) ; 
DirichletRandwert(ind) = 850*ones(size(ind,1),1);
ind = find( x(:,1) >= 6.691) ; 
DirichletRandwert(ind) = 850*ones(size(ind,1),1);
ind = find( x(:,1) >= 7.01) ; 
DirichletRandwert(ind) = 850*ones(size(ind,1),1);
ind = find( x(:,1) >= 7.51) ; 
DirichletRandwert(ind) = 850*ones(size(ind,1),1);
