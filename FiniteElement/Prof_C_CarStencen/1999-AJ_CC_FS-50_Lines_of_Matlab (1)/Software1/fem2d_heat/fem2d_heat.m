%FEM2D_HEAT   finite element method for two-dimensional heat equation.
%
%    FEM2D_HEAT solves heat equation 
%              d/dt u =  div(grad(u)) +  f  in Omega
%                   u = u_D                 on the Dirichlet boundary
%              d/dn u = g                   on the Neumann boundary
%    on a geometry described by triangles and presents the solution
%    graphically. 
% 
%    Therefore, FEM2D_HEAT assembles plane Courant finite elements and
%    calculates a discrete right-hand side as the coefficient vector of the
%    affine element approximation. Volume force and boundary data are given
%    as M-files <f.m>, <g.m>, and <u_d.m>. FEM2D_HEAT uses the reduced 
%    linear system of equations to calculate a discrete solution of the 
%    Laplace-problem. The resulting piecewise affine approximation will 
%    be graphically represented.
% 
%    FEM2D_heat loads the mesh data from data-files. The program reads the
%    triangular elements from the file <elements3.dat>.  The first column in
%    <elements3.dat> gives the number of each element. This is used for
%    clearness and is not neccesary for the numerical algorithm. The
%    following columns give the number of each node. Nodes of elements are
%    counted anti-clockwise.
% 
%    To adapt the program to a given heat equation the user has to specify
%    the data-files <coordinates.dat>, <elements3.dat>, <dirichlet.dat>, and
%    <neumann.dat> (optional) and the M-files <f.m>, <u_d.m>, and <g.m>
%    (optional). They have to be in the same directory as <fem2d_heat.m>.
%
%    Remark: This program is a supplement to the paper "Remarks around  
%    50 lines of Matlab: Short finite element implementation" by  
%    J. Alberty, C. Carstensen and S. A. Funken. The reader should 
%    consult that paper for more information.   
%
%
%    M-files you need to run FEM2D_HEAT
%       <stima3.m>, <f.m>, <u_d.m>, <show.m> and <g.m> (optional)
%
%    Data-files you need to run FEM2D_HEAT
%       <coordinates.dat>, <elements3.dat>, <dirichlet.dat>, and
%       <neumann.dat> (optional)

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <fem2d_heat.m> in $(HOME)/acf/fem2d_heat/
%    This program and corresponding data-files give Fig. 4d in 
%    "Remarks around 50 lines of Matlab: Short finite element 
%    implementation"

% Initialisation
load coordinates.dat; coordinates(:,1)=[];
load elements3.dat; elements3(:,1)=[];
eval('load neumann.dat; neumann(:,1) = [];','neumann=[];');
load dirichlet.dat; dirichlet(:,1) = [];
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse(size(coordinates,1),size(coordinates,1));
B = sparse(size(coordinates,1),size(coordinates,1));
T = 1; dt = 0.01; N = T/dt;
U = zeros(size(coordinates,1),N+1);

% Assembly
for j = 1:size(elements3,1)
  A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
      + stima3(coordinates(elements3(j,:),:));
end
for j = 1:size(elements3,1)
  B(elements3(j,:),elements3(j,:)) = B(elements3(j,:),elements3(j,:)) ...
      + det([1,1,1;coordinates(elements3(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
end

% Initial Condition
U(:,1) = zeros(size(coordinates,1),1);

% time steps
for n = 2:N+1
  b = sparse(size(coordinates,1),1);
  % Volume Forces
  for j = 1:size(elements3,1)
    b(elements3(j,:)) = b(elements3(j,:)) + ...
	det([1,1,1; coordinates(elements3(j,:),:)']) * ...
	dt*f(sum(coordinates(elements3(j,:),:))/3,n*dt)/6;
  end
  % Neumann conditions
  for j = 1 : size(neumann,1)
    b(neumann(j,:)) = b(neumann(j,:)) + ...
	norm(coordinates(neumann(j,1),:)-coordinates(neumann(j,2),:))* ...
	dt*g(sum(coordinates(neumann(j,:),:))/2,n*dt)/2; 
  end
  
  % previous timestep
  b = b + B * U(:,n-1);
  
  % Dirichlet conditions 
  u = sparse(size(coordinates,1),1);
  u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:),n*dt);
  b = b - (dt * A + B) * u;
  
  % Computation of the solution
  u(FreeNodes) = (dt*A(FreeNodes,FreeNodes)+ ...
      B(FreeNodes,FreeNodes))\b(FreeNodes); 
  U(:,n) = u;
end 

% graphic representation
show(elements3,[],coordinates,full(U(:,N+1)));
