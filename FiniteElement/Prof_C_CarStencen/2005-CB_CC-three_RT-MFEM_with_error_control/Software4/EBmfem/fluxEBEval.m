
function pEval=fluxEBEval(element,coordinate,u,nodes2edge,edge2element)

% ===========================================
% to compute p   by CB
% ===========================================
 weightedpoint = reshape(sum(reshape(coordinate(element',:),3,...
                  2*size(element,1))),size(element,1),2)/3;
 
 m=1;   
 pEval=zeros(3*size(element,1),2);   
 for j = 1:size(element,1)
   area(j) = det([1 1 1;coordinate(element(j,:),:)'])/2; 
   coord=coordinate(element(j,:),:) ; 
   %I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
   % check sign
   %tmp=signedge(find(j==edge2element(I,4)));
   I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
   signum=ones(1,3);   
   signum(find(j==edge2element(I,4)))=-1;
   c=coordinate(element(j,[2 3 1]),:)-coordinate(element(j,[3 1 2]),:);
   n=[norm(c(1,:)),norm(c(2,:)),norm(c(3,:))];
   
   % 2.way
   coord1=coordinate(element(j,:),:)'; 
   coordEval=coordinate(element(j,:),:)'*[4 1 1;1 4 1;1 1 4]/6;
   N=coordEval(:)*ones(1,3)-[coord1;coord1;coord1];
   PCC=N* 1/(2*area(j))*diag(signum)*diag(n)*u(I);
   PC=reshape(PCC,2,3);
   pEval(m:m+2,:)=[PC(1,:)',PC(2,:)'];
   m=m+3;  
   
   % 1.way
   % val = u(I)' * 1/(2*area(j))*diag(tmp)*diag(n);
   % ps1 = [4 1 1;1 4 1;1 1 4]/6 * coordinate(element(j,:),:) -  ones(3,1)*coord(1,:);
   % ps2 = [4 1 1;1 4 1;1 1 4]/6 * coordinate(element(j,:),:) -  ones(3,1)*coord(2,:);
   % ps3 = [4 1 1;1 4 1;1 1 4]/6 * coordinate(element(j,:),:)  -  ones(3,1)*coord(3,:);
	%ph = val(1,1)*ps1 + val(1,2)*ps2 + val(1,3)*ps3
   
end   


%% second way              % this gives the same result as the first way
  % coord1=coordinate(element(j,:),:)';   
  % Ans=coord1(:)*ones(1,3)
  % N=coord1(:)*ones(1,3)-[coord1;coord1;coord1];
   %PCC=((N* 1/(2*area(j))*diag(tmp)*diag(n))*u(I));
   %PC=reshape(PCC,2,3);
   %P(m:m+2,:)=[PC(1,:)',PC(2,:)'];