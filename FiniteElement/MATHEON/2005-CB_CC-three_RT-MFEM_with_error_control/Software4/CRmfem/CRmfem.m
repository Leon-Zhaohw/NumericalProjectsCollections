
% A. Program CRmfem for 2D Crouzeix-Thomas finite element method
%    
%
%    C.Bahriawati and C. Carstensen,08-08-03
%    File <CRmfem.m> 
%
%    M-files you need to run 
%       <stimaNC.m>, <edge.m>, <f.m>, <u_D.m>, <g.m> (optional)
%
%    Data-files you need to run 
%       <coordinate.dat>, <element.dat>, 
%       <Dirichlet.dat>, <Neumann.dat> (optional)
%
%    This program and corresponding data-files use Example 9.1 in
%    "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
%     MFEM with a Posteriori Error Control" by C.Bahriawati and C. Carstensen

% C.1. The main program 
load coordinate.dat;
load element.dat;
load Dirichlet.dat;
load Neumann.dat; 
[nodes2element,nodes2edge,noedges,edge2element,interioredge]=edge(element,coordinate);

% C.2. EBmfem
%function u = CRmfem(element,coordinate,Dirichlet,Neumann,...,
%                       nodes2element,nodes2edge,noedges,edge2element);
        
% Assemble matrix A
A=sparse(noedges,noedges); 
for j=1:size(element,1)
  I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
  A(I,I)=A(I,I)+StemaNC(coordinate(element(j,:),:));
end
% Volume force
b=sparse(noedges,1); 
for j=1:size(element,1)
  I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
  b(I,1)=b(I,1)+det([1 1 1;coordinate(element(j,:),:)'])*...
         f(sum(coordinate(element(j,:),:))/3)/6;
end
% Neumann conditions
if ~isempty(Neumann)
CN=coordinate(Neumann(:,2),:)-coordinate(Neumann(:,1),:);
for j = 1:size(Neumann,1)
  I=[nodes2edge(Neumann(j,1),Neumann(j,2))];
  edgeNeumann=nodes2edge(Neumann(j,1),Neumann(j,2));
  b(edgeNeumann,1)=b(edgeNeumann,1)+ norm(CN(j,:))*...
                 g(sum(coordinate(Neumann(j,:),:))/2,CN(j,:)*[0,-1;1,0]/norm(CN(j,:)));             
 end
end
% Dirichlet conditions
tmp=zeros(noedges,1);
tmp(diag(nodes2edge(Dirichlet(:,1),Dirichlet(:,2))))=...
    ones(size(diag(nodes2edge(Dirichlet(:,1),Dirichlet(:,2))),1),1);
FreeNodes=find(~tmp);
x=zeros(noedges,1);
x(diag(nodes2edge(Dirichlet(:,1),Dirichlet(:,2))))=...
  u_D((coordinate(Dirichlet(:,1),:)+coordinate(Dirichlet(:,2),:))/2);
b=b-A*x;
x(FreeNodes)=A(FreeNodes,FreeNodes)\b(FreeNodes);

u=x;

% Compute the fluxes
ph=ph_OnRTElement(element,coordinate,nodes2edge,noedges,...
                  edge2element,u);
% ==========================================
% Compute uh on L_0(T) ( for CR-ThomasElemnt)
% ===========================================            
calculated=zeros(size(element,1),1);
IntPh=IntPhOmega(element,coordinate,u,nodes2edge,edge2element);
[uhD,calculated]=uhDir(coordinate,element,Dirichlet,nodes2edge,...
                      nodes2element,edge2element,ph,calculated,u);           
uh=uhD; 
[uh,calculated]=uhN(element,coordinate,nodes2edge,edge2element,...
               noedges,nodes2element,ph,uh,IntPh,calculated); 
figure(1)
ShowDisplacement(element,coordinate,uh);  % u and uh
figure(2)
ShowFlux(element,coordinate,ph);
           

           
