%FEM2D_NONLINEAR  finite element method for two-dimensional nonlinear equation.
%
%    FEM2D_NONLINEAR solves the two-dimensional nonlinear equation 
%     K(u) + f = 0 in Omega
%            u = u_D  on the Dirichlet boundary
%       d/dn u = g   on the Neumann boundary
%    on a geometry described by triangles and presents a solution
%    graphically. (The solutions are in general not unique.)
% 
%    Therefore, FEM2D_NONLINEAR solves the weak problem J(u,v) + <f,v> = 0
%    with Newton-Raphson's method. Starting with some initial value U in
%    each iteration step FEM2D_NONLINEAR solves the discretized version of
%    the linearized Problem
%
%           DJ(U) * W = J(U) + f
%
%    for W and then updates U as U = U - W.
%
%    DJ(U) is the linearization of J at U. J resp. DJ are given as M-files
%    <localj.m> resp. <localdj.m>. Volume force and boundary data are given 
%    as M-files <f.m>, <g.m>, and <u_d.m>. FEM2D_NONLINEAR uses the reduced 
%    linear system of equations to calculate a discrete solution of the 
%    linearized problem. The number of iteration steps is limited to 50. The 
%    resulting piecewise affine approximation will be graphically represented.
% 
%    FEM2D_NONLINEAR loads the mesh data from data-files. The program reads
%    the triangular elements from the file <elements3.dat>.  The first column
%    in <elements3.dat> gives the number of each element. This is used for
%    clearness and is not neccesary for the numerical algorithm. The
%    following columns give the number of each node. Nodes of elements are
%    counted anti-clockwise.
% 
%    To adapt the program to a given nonlinear equation the user has to
%    specify the data-files <coordinates.dat>, <elements3.dat>,
%    <dirichlet.dat>, and <neumann.dat> (optional) and the M-files <f.m>,
%    <u_d.m>, and <g.m> (optional). They have to be in the same directory as
%    <fem2d_nonlinear.m>.
%
%    Remark: This program is a supplement to the paper "Remarks around  
%    50 lines of Matlab: Short finite element implementation" by  
%    J. Alberty, C. Carstensen and S. A. Funken. The reader should 
%    consult that paper for more information.   
%
%
%    M-files you need to run FEM2D_NONLINEAR
%       <localj.m>, <localdj.m>, <f.m>, <u_d.m>, <show.m> and <g.m> (optional)
%
%    Data-files you need to run FEM2D_NONLINEAR
%       <coordinates.dat>, <elements3.dat>,
%       <dirichlet.dat>, and <neumann.dat> (optional)

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <fem2d_nonlinear.m> in $(HOME)/acf/fem2d_nonlinear/
%    This program and corresponding data-files give Fig. 5a in 
%    "Remarks around 50 lines of Matlab: Short finite element 
%    implementation"

%
% Initialisation
load coordinates.dat; coordinates(:,1)=[];
load elements3.dat; elements3(:,1)=[];
eval('load neumann.dat; neumann(:,1) = [];','neumann=[];');
load dirichlet.dat; dirichlet(:,1) = [];
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));

% Initial value
U = -ones(size(coordinates,1),1);
U(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));

% Newton-Raphson iteration
for i=1:50
  
  % Assembly of DJ(U)
  A = sparse(size(coordinates,1),size(coordinates,1));
  for j = 1:size(elements3,1)
    A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
	+ localdj(coordinates(elements3(j,:),:),U(elements3(j,:)));
  end
  
  % Assembly of J(U)
  b = sparse(size(coordinates,1),1);
  for j = 1:size(elements3,1)
    b(elements3(j,:)) = b(elements3(j,:)) ...
   	+ localj(coordinates(elements3(j,:),:),U(elements3(j,:)));
  end
  
  % Volume Forces
  for j = 1:size(elements3,1)
    b(elements3(j,:)) = b(elements3(j,:)) + ...
	det([1 1 1; coordinates(elements3(j,:),:)']) * ...
	f(sum(coordinates(elements3(j,:),:))/3)/6;
  end
  
  % Neumann conditions
  for j = 1 : size(neumann,1)
    b(neumann(j,:))=b(neumann(j,:)) - norm(coordinates(neumann(j,1),:)- ...
	coordinates(neumann(j,2),:))*g(sum(coordinates(neumann(j,:),:))/2)/2;
  end
  
  % Dirichlet conditions
  W = zeros(size(coordinates,1),1);
  W(unique(dirichlet)) = 0;
  
  % Solving one Newton step
  W(FreeNodes) = A(FreeNodes,FreeNodes)\b(FreeNodes);
  U = U - W;
  if norm(W) < 10^(-10)
    break
  end
end

% graphic representation
show(elements3,[],coordinates,full(U));
