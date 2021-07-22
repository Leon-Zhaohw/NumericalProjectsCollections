function [AvE,Eps3,Eps4,AvS,Sigma3,Sigma4] = ...
                     avmatrix(coordinates,elements3,elements4,u,lambda,mu);
%AVMATRIX  preprocessing for estimator and visualisation.
%   [AVE,EPS3,EPS4,AVS,SIGMA3,SIGMA4] = AVMATRIX(COORDINATES,ELEMENTS3,
%   ELEMENTS4,U,LAMBDA,MU) determines variables used as input for the
%   functions <aposteriori.m> and <show.m>.
%
%
%   See also FEM_LAME2D, AVMATRIX, and SHOW.

%    J. Alberty, C. Carstensen and S. A. Funken  07-03-00
%    File <avmatrix.m> in $(HOME)/acfk/fem_lame2d/cooks/ and
%                         $(HOME)/acfk/fem_lame2d/lshape_p1/ and
%                         $(HOME)/acfk/fem_lame2d/lshape_q1/ and
%                         $(HOME)/acfk/fem_lame2d/hole/

Eps3 = zeros(size(elements3,1),4);
Sigma3 = zeros(size(elements3,1),4);
Eps4 = zeros(size(elements4,1),4);
Sigma4 = zeros(size(elements4,1),4);
AreaOmega = zeros(size(coordinates,1),1);
AvS = zeros(size(coordinates,1),4);
AvE = zeros(size(coordinates,1),4);
for j = 1:size(elements3,1)
  area3 = det([1,1,1;coordinates(elements3(j,:),:)'])/2;
  AreaOmega(elements3(j,:)) = AreaOmega(elements3(j,:)) +area3;
  PhiGrad = [1,1,1;coordinates(elements3(j,:),:)']\[zeros(1,2);eye(2)];
  U_Grad = u([1;1]*2*elements3(j,:)-[1;0]*[1,1,1])*PhiGrad;
  Eps3(j,:) = reshape((U_Grad+U_Grad')/2,1,4);
  Sigma3(j,:) = reshape(lambda*trace(U_Grad)*eye(2) ...
      +2*mu*(U_Grad+U_Grad')/2,1,4);
  AvE(elements3(j,:),:) = AvE(elements3(j,:),:) +area3*[1;1;1]*Eps3(j,:);
  AvS(elements3(j,:),:) = AvS(elements3(j,:),:) +area3*[1;1;1]*Sigma3(j,:);
end;
for j = 1:size(elements4,1)
  area4 = det([1,1,1;coordinates(elements4(j,1:3),:)']);
  AreaOmega(elements4(j,:),:) = AreaOmega(elements4(j,:)) +area4;
  Vertices = coordinates(elements4(j,:),:);
  U_Grad = inv([Vertices(2,:)-Vertices(1,:);Vertices(4,:)-...
	Vertices(1,:)])*[-1,1,1,-1;-1,-1,1,1]*...
      [u(2*elements4(j,:)-1),u(2*elements4(j,:))]/2;
  Eps4(j,:) = reshape((U_Grad+U_Grad')/2,1,4);
  Sigma4(j,:) = reshape(lambda*trace(U_Grad)*eye(2) ...
      +2*mu*(U_Grad+U_Grad')/2,1,4);
  AvE(elements4(j,:),:) = AvE(elements4(j,:),:) +area4*[1;1;1;1]*Eps4(j,:);
  AvS(elements4(j,:),:) = AvS(elements4(j,:),:) +area4*[1;1;1;1]*Sigma4(j,:);
end
AvE = AvE./(AreaOmega*[1,1,1,1]);
AvS = AvS./(AreaOmega*[1,1,1,1]);
