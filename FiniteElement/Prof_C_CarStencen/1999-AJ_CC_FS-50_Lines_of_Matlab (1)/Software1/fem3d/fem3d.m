%FEM3D  three-dimensional finite element method for Laplacian.
%
%    FEM3D solves Laplace's equation 
%      - div(grad(u)) = f in Omega
%                   u = u_D on the Dirichlet boundary
%              d/dn u = g   on the Neumann boundary
%    on a geometry described by tetraeder and presents the solution 
%    graphically.
% 
%    Therefore, FEM3D assembles picewise linear finite elements and 
%    calculates a discrete right-hand side as the coefficient vector  
%    of the affine element approximation. Volume force and boundary 
%    data are given as M-files <f.m>, <g.m>, <u_d.m>. FEM3D uses 
%    the reduced linear system of equations to calculate a discrete 
%    solution of the Laplace-problem. The resulting piecewise affine
%    approximation will be graphically represented. 
% 
%    FEM3D loads the mesh data from data-files. The program reads 
%    the tetraeder elements from the files <elements3.dat>. 
%    The first column in <elements3.dat> gives the 
%    number of each element. This is used for clearness and 
%    is not neccesary for the numerical algorithm. The following 
%    columns give the number of each node. 
% 
%    To adapt the program to a given Laplace equation the user has to 
%    specify the data-files <coordinates.dat>, <elements3.dat>, 
%    <dirichlet.dat>, <neumann.dat> (optional) 
%    and the M-files <f.m>, <g.m>, <u_d.m>. They have to be in the 
%    same directory as <fem3d.m>. 
%
%    Remark: This program is a supplement to the paper "Remarks around  
%    50 lines of Matlab: Short finite element implementation" by  
%    J. Alberty, C. Carstensen and S. A. Funken. The reader should 
%    consult that paper for more information.   
%
%
%    M-files you need to run FEM3D
%       <stima3.m>, <f.m>, <u_d>, <showsurface.m> and <g.m> (optional)
%
%    Data-files you need to run FEM3D
%       <coordinates.dat>, <elements3.dat>, 
%       <dirichlet.dat>, <neumann.dat> (optional)

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <fem3d.m> in $(HOME)/acf/fem3d/
%    This program and corresponding data-files give Fig. 6 in 
%    "Remarks around 50 lines of Matlab: Short finite element 
%    implementation"

% Initialisation
load coordinates.dat; coordinates(:,1)=[];
load elements3.dat; elements3(:,1)=[];
eval('load neumann.dat; neumann(:,1) = [];','neumann=[];');
load dirichlet.dat; dirichlet(:,1) = [];
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse(size(coordinates,1),size(coordinates,1));
b = sparse(size(coordinates,1),1);

% Assembly
for j = 1:size(elements3,1)
  A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
      + stima3(coordinates(elements3(j,:),:));
end

% Volume Forces
for j = 1:size(elements3,1)
  b(elements3(j,:)) = b(elements3(j,:)) + ...
      det([1,1,1,1;coordinates(elements3(j,:),:)']) ... 
      * f(sum(coordinates(elements3(j,:),:))/4) / 24;
end

% Neumann conditions
for j = 1 : size(neumann,1)
  b(neumann(j,:)) = b(neumann(j,:)) + ...
      norm(cross(coordinates(neumann(j,3),:)-coordinates(neumann(j,1),:), ...
      coordinates(neumann(j,2),:)-coordinates(neumann(j,1),:))) ...
      * g(sum(coordinates(neumann(j,:),:))/3)/6;
end

% Dirichlet conditions 
u = sparse(size(coordinates,1),1);
u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));
b = b - A * u;

% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

% Graphic representation
showsurface([dirichlet;neumann],coordinates,full(u));
