% A. Program EBMFEM for 2D Raviart-Thomas mixed finite element method
%    based on the edge-oriented basis function
%
%    C.Bahriawati and C. Carstensen,08-08-03
%    File <EBmfem.m> 
%
%    M-files you need to run 
%       <stimaB.m>, <edge.m>, <f.m>, <u_D.m>, <g.m> (optional)
%
%    Data-files you need to run 
%       <coordinate.dat>, <element.dat>, 
%       <Dirichlet.dat>, <Neumann.dat> (optional)
%
%    This program and corresponding data-files use Example 9.1 in
%    "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
%     MFEM with a Posteriori Error Control" by C.Bahriawati and C. Carstensen

% A.1. The main program 
load coordinate.dat;
load element.dat;
load dirichlet.dat;
load Neumann.dat; 
[nodes2element,nodes2edge,noedges,edge2element,interioredge]=edge(element,coordinate);

% A.2. EBmfem
%function u=EBmfem(element,coordinate,dirichlet,Neumann,nodes2element,...
%         nodes2edge,noedges,edge2element);

% Assemble matrices B and C    
 B=sparse(noedges, noedges); 
 C=sparse(noedges,size(element,1));
for j = 1:size(element,1)
   coord=coordinate(element(j,:),:)';
   I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
   signum=ones(1,3);   
   signum(find(j==edge2element(I,4)))=-1;
   B(I,I)= B(I,I)+diag(signum)*stimaB(coord)*diag(signum); 
   n=coord(:,[3,1,2])-coord(:,[2,3,1]);
   C(I,j) = diag(signum)*[norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]';
end
% Global stiffness matrix A
  A = sparse(noedges+size(element,1), noedges+size(element,1));
  A = [B ,         C,      ; 
       C', sparse(size(C,2),size(C,2))];    
% Volume force 
b = sparse(noedges+size(element ,1),1);
for j = 1:size(element ,1)    
  b(noedges+j)= -det([1,1,1; coordinate(element(j,:),:)']) * ...
                f(sum(coordinate(element(j,:),:))/3)/6;
end
% Dirichlet conditions
for k = 1:size(dirichlet,1)
  b(nodes2edge(dirichlet(k,1),dirichlet(k,2)))= norm(coordinate(dirichlet(k,1),:)-...
       coordinate(dirichlet(k,2),:))*u_D(sum(coordinate(dirichlet(k,:),:))/2);
end   
% Neumann conditions
if ~isempty(Neumann)
tmp=zeros(noedges+size(element,1),1);
tmp(diag(nodes2edge(Neumann(:,1),Neumann(:,2))))=...
   ones(size(diag(nodes2edge(Neumann(:,1),Neumann(:,2))),1),1);
FreeEdge=find(~tmp);
x=zeros(noedges+size(element,1),1);
CN=coordinate(Neumann(:,2),:)-coordinate(Neumann(:,1),:);
for j=1:size(Neumann,1)
 x(nodes2edge(Neumann(j,1),Neumann(j,2)))=...
 g(sum(coordinate(Neumann(j,:),:))/2,CN(j,:)*[0,-1;1,0]/norm(CN(j,:)));
end
b=b-A*x;
x(FreeEdge)=A(FreeEdge,FreeEdge)\b(FreeEdge);   
else
    x = A\b;
end     
figure(1)
ShowDisplacement(element,coordinate,x);
p=fluxEB(element,coordinate,x,noedges,nodes2edge,edge2element);
figure(2)
ShowFlux(element,coordinate,p);
pEval=fluxEBEval(element,coordinate,x,nodes2edge,edge2element);
eta_T = Aposteriori(element,coordinate,dirichlet,Neumann,x,pEval)             





