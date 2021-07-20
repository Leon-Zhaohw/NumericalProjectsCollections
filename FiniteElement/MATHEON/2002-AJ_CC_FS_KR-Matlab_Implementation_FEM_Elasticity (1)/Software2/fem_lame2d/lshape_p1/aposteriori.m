function estimate=aposteriori(coordinates,elements3,elements4, ...
                                AvE,Eps3,Eps4,AvS,Sigma3,Sigma4,u,lambda,mu);
%APOSTERIORI  a posteriori error estimator.
%   ETA = APOSTERIORI(COORDINATES,ELEMENTS3,ELEMENTS4,AVE,EPS3,EPS4,AVS,
%   SIGMA3,SIGMA4,U,LAMBDA,MU) returns the estimated error of the numerical
%   solution U for the mesh defined by COORDINATES, ELEMENTS3, and
%   ELEMENTS4. LAMBDA and MU are the Lame constants, the variables AVE,
%   EPS3, EPS4, AVS, SIGMA3, and SIGMA4 are previously determined by the
%   function <avmatrix.m>.
%
%
%   See also FEM_LAME2D and AVMATRIX.

%    J. Alberty, C. Carstensen and S. A. Funken  07-03-00
%    File <aposteriori.m> in $(HOME)/acfk/fem_lame2d/cooks/ and
%                            $(HOME)/acfk/fem_lame2d/lshape_p1/ and
%                            $(HOME)/acfk/fem_lame2d/lshape_q1/ and
%                            $(HOME)/acfk/fem_lame2d/hole/

eta3 = zeros(size(elements3,1),1);
eta4 = zeros(size(elements4,1),1);
e3 = zeros(size(elements3,1),1);
e4 = zeros(size(elements4,1),1);
for j=1:size(elements3,1)
  Area3=det([1,1,1;coordinates(elements3(j,:),:)'])/2;
  for i=1:4 
    eta3(j) = eta3(j) + Area3 * (Sigma3(j,i)*Eps3(j,i) ...
	+ AvS(elements3(j,:),i)'*[2,1,1;1,2,1;1,1,2]*AvE(elements3(j,:),i)/12 ...
	- AvS(elements3(j,:),i)'*[1;1;1]*Eps3(j,i)/3 ...
	- AvE(elements3(j,:),i)'*[1;1;1]*Sigma3(j,i)/3);
    e3(j) = e3(j) + Area3 * ...
	AvS(elements3(j,:),i)'*[2,1,1;1,2,1;1,1,2]*AvE(elements3(j,:),i)/12;
  end
end
for j=1:size(elements4,1)
  Area4=det([1,1,1;coordinates(elements4(j,1:3),:)']);
  for i=1:4
    eta4(j) = eta4(j) + Area4 * (Sigma4(j,i)*Eps4(j,i) ...
	+ AvS(elements4(j,:),i)'*[4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4]* ...
	AvE(elements4(j,:),i)/36 ...
	- AvS(elements4(j,:),i)'*[1;1;1;1]*Eps4(j,i)/4 ...
	- AvE(elements4(j,:),i)'*[1;1;1;1]*Sigma4(j,i)/4);
    e4(j) = e4(j) + Area4 * ...
	AvS(elements4(j,:),i)'*[4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4]* ...
	AvE(elements4(j,:),i)/36;
  end
end
estimate = sqrt(sum(eta3)+sum(eta4))/ sqrt(sum(e3)+sum(e4));
