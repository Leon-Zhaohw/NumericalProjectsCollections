function VolumeForce = f(x);
%F   Volume force in considered domain.
%   Y = F(X) returns values of forces at N discrete points in the considered
%   domain. This input data has to be chosen by the user. X has dimension N
%   x 2 and Y has dimension N x 1.
%
%
%   See also FEM2D, U_D, and G.

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <f.m> in $(HOME)/acf/fem2d/
%    This volume force is used to compute Fig. 3 in 
%    "Remarks around 50 lines of Matlab: Short finite element 
%    implementation".

VolumeForce = ones(size(x,1),1);
