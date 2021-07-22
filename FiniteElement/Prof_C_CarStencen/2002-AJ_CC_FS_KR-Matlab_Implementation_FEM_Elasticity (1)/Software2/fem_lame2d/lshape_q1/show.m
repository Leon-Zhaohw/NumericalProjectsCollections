function show(elements3,elements4,coordinates,AvS,u,lambda,mu)
%SHOW  Plots two-dimensional solution
%    SHOW(ELEMENTS3,ELEMENTS4,COORDINATES,AVS,U) plots the
%    strained mesh and visualizes the stresses in grey tones.
%
%    The variable AVS is previously determined by the function <avmatrix.m>.
%
%
%   See also FEM_LAME2D and AVMATRIX.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <show.m> in $(HOME)/acfk/fem_lame2d/cooks/ and
%                     $(HOME)/acfk/fem_lame2d/lshape_p1/ and
%                     $(HOME)/acfk/fem_lame2d/lshape_q1/ and
%                     $(HOME)/acfk/fem_lame2d/hole/

for i=1:size(coordinates,1)
 AvC(i)=(mu/(24*(mu+lambda)^2)+1/(8*mu))*(AvS(i,1)+...
        AvS(i,4))^2+1/(2*mu)*(AvS(2)^2-AvS(1)*AvS(4));
end
factor=20;
colormap(1-gray)
trisurf(elements3,factor*u(1:2:size(u,1))+coordinates(:,1), ...
    factor*u(2:2:size(u,1))+coordinates(:,2), ...
    zeros(size(coordinates,1),1), AvC, 'facecolor','interp');
hold on
trisurf(elements4,factor*u(1:2:size(u,1))+coordinates(:,1), ...
    factor*u(2:2:size(u,1))+coordinates(:,2), ...
    zeros(size(coordinates,1),1), AvC, 'facecolor','interp');
view(0,90)
hold off
colorbar('vert')
  



