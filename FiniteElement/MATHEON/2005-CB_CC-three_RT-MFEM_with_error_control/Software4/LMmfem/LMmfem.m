% B. Program LMMFEM for 2D Raviart-Thomas mixed finite element method
%    based on the Lagrange Multiplier Technique
%
%    C.Bahriawati and C. Carstensen,08-08-03
%    File <LMmfem.m> 
%
%    M-files you need to run 
%       <edge.m>, <f.m>, <u_D.m>, <g.m> (optional)
%
%    Data-files you need to run 
%       <coordinate.dat>, <element.dat>, 
%       <Dirichlet.dat>, <Neumann.dat> (optional)
%
%    This program and corresponding data-files use Example 9.1 in
%    "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
%     MFEM with a Posteriori Error Control" by C.Bahriawati and C. Carstensen

% B.1. The main program 
load coordinate.dat;
load element.dat;
load dirichlet.dat;
load Neumann.dat; 
[nodes2element,nodes2edge,noedges,edge2element,interioredge]=edge(element,coordinate);

% B.2. LMmfem
%function u=LMmfem(element,coordinate,dirichlet,Neumann,nodes2element,interioredge);
MidPoint=reshape(sum(reshape(coordinate(element',:),3,...
         2*size(element,1))),size(element,1),2)/3;
%Assemble matrices B and C
B=sparse(3*size(element,1),3*size(element,1));
C=sparse(3*size(element,1),size(element,1));
for j=1:size(element,1)
  s=sum(sum((coordinate(element(j,[2 3 1]),:)- ...
                coordinate(element (j,[1 2 3]),:)).^2));
  B(3*(j-1)+[1 2 3],3*(j-1)+[1 2 3])=det([1 1 1;coordinate(element (j,:),:)'])*...
                                     diag([1 1 s/36])/2;
  C(3*(j-1)+[1,2,3],j)=[0;0;det([1 1 1;coordinate(element (j,:),:)'])];                               
end
% Assemble matrix D                                                                               
D=sparse(3*size(element ,1),size(interioredge,1));
EV=coordinate(interioredge(:,2),:)-coordinate(interioredge(:,1),:);                       
for k=1:size(interioredge,1)                                                            
  h=(coordinate(interioredge(k,[1,1]),:)-MidPoint(interioredge(k,3:4),:))*...
                                                          [EV(k,2);-EV(k,1)];                                                                                                        
  D([3*(interioredge(k,3)-1)+[1,2,3],3*(interioredge(k,4)-1)+[1,2,3]],k)=...           
                            [-EV(k,2);EV(k,1);-h(1); EV(k,2);-EV(k,1); h(2)];           
end
% Global stiffness matrix A                                                             
 A = sparse(4*size(element,1)+size(interioredge,1),...
            4*size(element,1)+size(interioredge,1));                            
 A = [B ,         C,              D           ;
      C',sparse(size(C,2),size(C,2)+size(D,2));
      D',sparse(size(D,2),size(C,2)+size(D,2))];
% Volume force
b=sparse(4*size(element,1)+size(interioredge,1),1);
for j=1:size(element,1)
  b(3*size(element,1)+j)=-det([1,1,1;coordinate(element(j,:),:)'])*...
                          f(sum(coordinate(element(j,:),:))/3)/6; 
end
% Dirichlet conditions
EV=coordinate(dirichlet(:,2),:)-coordinate(dirichlet(:,1),:);
for k=1:size(dirichlet,1)   
  ElementDir=nodes2element(dirichlet(k,1),dirichlet(k,2));
  h=(coordinate(dirichlet(k,1),:)-MidPoint(ElementDir,:))*...
     [EV(k,2);-EV(k,1)]; 
  b(3*nodes2element(dirichlet(k,1),dirichlet(k,2))-[2 1 0])=...
  b(3*nodes2element(dirichlet(k,1),dirichlet(k,2))-[2 1 0]) + ...
    u_D(sum(coordinate(dirichlet(k,:),:))/2)*[EV(k,2);-EV(k,1);h]; 
end
% Neumann conditions
if ~isempty(Neumann)
   F=sparse(3*size(element,1),size(Neumann,1));
   CN=coordinate(Neumann(:,2),:) - coordinate(Neumann(:,1),:);
   for k=1:size(Neumann,1)  
    h=(coordinate(Neumann(k,1),:)-...
       MidPoint(nodes2element(Neumann(k,1),Neumann(k,2)),:))*[CN(k,2);-CN(k,1)];    
       F([3*(nodes2element(Neumann(k,1),Neumann(k,2))-1)+[1,2,3]],k)=...
                                                           [CN(k,2);-CN(k,1);h];                       
 end
 F=[F;sparse(size(A,1)-size(F,1),size(F,2))];
 A=[A,F;F',sparse(size(F,2),size(F,2))];
 % Right-hand side  
 b=[b;sparse(size(Neumann,1),1)];
 for j=1:size(Neumann,1)
   b(4*size(element,1)+size(interioredge,1)+j)=...
   b(4*size(element,1)+size(interioredge,1)+j) + ...
   norm(CN(j,:))*g(sum(coordinate(Neumann(j,:),:))/2,CN(j,:)*[0,-1;1,0]/norm(CN(j,:)));
 end
end
% Computation of the solution
x=A\b;

figure(1)
uinit=3*size(element,1)+1:4*size(element,1);
ShowDisplacement(element,coordinate,x(uinit));
figure(2)
p=fluxLM(element,coordinate,x);
ShowFlux(element,coordinate,p);
pEval=fluxLMEval(element,coordinate,x);
eta_T = Aposteriori(element,coordinate,dirichlet,Neumann,x,pEval)







