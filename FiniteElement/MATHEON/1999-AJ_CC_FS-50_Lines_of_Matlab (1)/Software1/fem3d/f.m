function VolumeForce = f(x);
%F   Volume force in considered domain.
%   Y = F(X) returns values of forces at N discrete points in the considered
%   domain. This input data has to be chosen by the user. X has dimension N
%   x 3 and Y has dimension N x 1.
%
%
%   See also FEM3D, U_D, and G.

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <f.m> in $(HOME)/acf/fem3d/
%    This volume force is used to compute Fig. 6 in 
%    "Remarks around "50 lines of Matlab: Short finite element 
%    implementation".

VolumeForce = 0.001*x(:,1)^3;
