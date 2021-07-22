% C.Bahriawati
% Compute Div \psi on T
% ==================================


function [ps,val]=DivPsi(coordinate,element,nodes2edge,edge2element,N)

coord=coordinate(element(N,:),:);
Area = det([1 1 1;coordinate(element(N,:),:)'])/2; 
L=diag(nodes2edge(element(N,[2 3 1]),element(N,[3 1 2])));
     % Check sign
[r,s]=find(N==edge2element(L,4));
tmp=signedge(r);
c=coordinate(element(N,[2 3 1]),:)-coordinate(element(N,[3 1 2]),:);
n=[norm(c(1,[1 2])),norm(c(2,[1 2])),norm(c(3,[1 2]))];
     % ==========
val = 1/(2*Area)*diag(tmp)*diag(n);
ps=repmat([4 1 1;1 4 1;1 1 4]/6 * coordinate(element(N,:),:),3,1)-...
   [ones(3,1)*coord(1,:);ones(3,1)*coord(2,:);ones(3,1)*coord(3,:)];    
    
    
     
     % 1.way 
    %ps1 = [4 1 1;1 4 1;1 1 4]/6 * coordinate(element(N,:),:) -  ones(3,1)*coord(1,:);
    %ps2 = [4 1 1;1 4 1;1 1 4]/6 * coordinate(element(N,:),:) -  ones(3,1)*coord(2,:);
    %ps3 = [4 1 1;1 4 1;1 1 4]/6 * coordinate(element(N,:),:) -  ones(3,1)*coord(3,:);
    %ps=[ps1;ps2;ps3]
     
     % 2.way 
   %coord1=coordinate(element(N,:),:)';
   %coordEval1=coordinate(element(N,:),:)'*[4 1 1;1 4 1;1 1 4]/6;
   %N=coordEval1(:)*ones(1,3)-[coord1;coord1;coord1];
   %PCC=N* 1/(2*area(ElDir))*diag(tmp)*diag(n)*uI;
   %PC=reshape(PCC,2,3)
   %pEval=[PC(1,:)',PC(2,:)'] - (fbar.*ps1+fbar.*ps2+fbar.*ps3)