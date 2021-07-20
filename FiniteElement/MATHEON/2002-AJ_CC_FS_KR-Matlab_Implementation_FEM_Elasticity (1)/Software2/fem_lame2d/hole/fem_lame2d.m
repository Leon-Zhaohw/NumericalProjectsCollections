%FEM_LAME2D  two-dimensional finite element method for the Lame-Problem
%
%    FEM_LAME2D solves the Lame Problem
%      (lambda+mu)(grad div u) +mu div grad u +f = 0  in Omega
%                                            M u = w  on the Dirichlet boundary
%           (lambda tr eps(u) Id +2 mu eps(u)) n = g  on the Neumann boundary
%    on a geometry described by triangles and parallelograms and presents
%    the solution graphically.
%
%    The parameters E and NU, defined in the first line, are problem
%    dependend and have to be chosen by the user.
%
%    Remark: This program is a supplement to the paper
%    "Matlab-Implementation of the Finite Element Method in Elasticity" by  
%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose. The reader should 
%    consult that paper for more information.   
%
%
%    M-files you need to run FEM_LAME2D
%       <stima3.m>, <stima4.m>, <f.m>, <u_d.m>, <avmatrix.m>, <show.m>,
%       <aposteriori.m>, and <g.m> (optional)
%
%    Data-files you need to run FEM_LAME2D
%       <coordinates.dat>, <elements.dat>, 
%       <dirichlet.dat>, <neumann.dat> (optional)

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <fem_lame2d.m> in $(HOME)/acfk/fem_lame2d/cooks/ and
%                           $(HOME)/acfk/fem_lame2d/hole/
%    This program and corresponding data-files give Fig. 6 resp. Fig. 7 in
%    "Matlab-Implementation of the Finite Element Method in Elasticity"

% Initialisation
E = 2900; nu = 0.4;
mu = E/(2*(1+nu)); lambda = E*nu/((1+nu)*(1-2*nu));
load coordinates.dat;
eval('load elements3.dat;','elements3 = [];');
eval('load elements4.dat;','elements4 = [];');
eval('load neumann.dat;','neumann = [];');
load dirichlet.dat;
A = sparse(2*size(coordinates,1),2*size(coordinates,1)); 
b = zeros(2*size(coordinates,1),1);

% Assembly
for j = 1:size(elements3,1)
  I = 2*elements3(j,[1,1,2,2,3,3]) -[1,0,1,0,1,0]; 
  A(I,I) = A(I,I) +stima3(coordinates(elements3(j,:),:),lambda,mu);   
end
for j = 1:size(elements4,1)
  I = 2*elements4(j,[1,1,2,2,3,3,4,4]) -[1,0,1,0,1,0,1,0];
  A(I,I) = A(I,I) +stima4(coordinates(elements4(j,:),:),lambda,mu);   
end

% Volume forces
for j = 1:size(elements3,1)
  I = 2*elements3(j,[1,1,2,2,3,3]) -[1,0,1,0,1,0];
  fs = f(sum(coordinates(elements3(j,:),:))/3)';
  b(I) = b(I) +det([1,1,1;coordinates(elements3(j,:),:)'])*[fs;fs;fs]/6;
end
for j = 1:size(elements4,1)
  I = 2*elements4(j,[1,1,2,2,3,3,4,4]) -[1,0,1,0,1,0,1,0];
  fs = f(sum(coordinates(elements4(j,:),:))/4)';
  b(I) = b(I) +det([1,1,1;coordinates(elements4(j,1:3),:)'])*[fs;fs;fs;fs]/4;
end

% Neumann conditions
if ~isempty(neumann)
  n = (coordinates(neumann(:,2),:) -coordinates(neumann(:,1),:))*[0,-1;1,0];
  for j = 1:size(neumann,1);
    I = 2*neumann(j,[1,1,2,2]) -[1,0,1,0];
    gm = g(sum(coordinates(neumann(j,:),:))/2, n(j,:)/norm(n(j,:)))';
    b(I) = b(I) +norm(n(j,:))*[gm;gm]/2;
  end
end

% Dirichlet conditions
DirichletNodes = unique(dirichlet);
[W,M] = u_d(coordinates(DirichletNodes,:));
B = sparse(size(W,1),2*size(coordinates,1));
for k = 0:1
  for l = 0:1
    B(1+l:2:size(M,1),2*DirichletNodes-1+k) = diag(M(1+l:2:size(M,1),1+k));
  end
end
mask = find(sum(abs(B)'));
A = [A, B(mask,:)'; B(mask,:), sparse(length(mask),length(mask))];
b = [b;W(mask,:)];

% Calculating the solution
x = A \ b;
u = x(1:2*size(coordinates,1));

% Representation of the solution
[AvE,Eps3,Eps4,AvS,Sigma3,Sigma4] = ...
    avmatrix(coordinates,elements3,elements4,u,lambda,mu);
show(elements3,elements4,coordinates,AvS,u,lambda,mu);
estimate = aposteriori(coordinates,elements3,elements4,AvE,Eps3,Eps4, ...
    AvS,Sigma3,Sigma4,u,lambda,mu)
