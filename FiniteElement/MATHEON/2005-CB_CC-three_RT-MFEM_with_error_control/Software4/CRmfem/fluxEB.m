
function p=fluxEB(element,coordinate,u,noedges,nodes2edge,edge2element)
%function p=fluxEB(element,coordinate,u,noedges,nodes2edge,edge2element);

% %%%%%%%%%%%%%%%%%%%
 % to compute p
% %%%%%%%%%%%%%%%%%%%
 weightedpoint = reshape(sum(reshape(coordinate(element',:),3,...
                  2*size(element,1))),size(element,1),2)/3;
 
 %m=1;
 p=zeros(3*size(element,1),2);   
 P=zeros(3*size(element,1),2);   
 
 for j = 1:size(element,1)
   area(j) = det([1 1 1;coordinate(element(j,:),:)'])/2; 
   coord=coordinate(element(j,:),:);
   I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
   % check sign
   %r=find(j==edge2element(I,4));
   %tmp=signedge(find(j==edge2element(I,4)));
   signum=ones(1,3);  % Edited version 7/4/3: 1 pm
   signum(find(j==edge2element(I,4)))=-1;
   c=coordinate(element(j,[2 3 1]),:)-coordinate(element(j,[3 1 2]),:);
   n=[norm(c(1,:)),norm(c(2,:)),norm(c(3,:))];
   %% first way
   p1 = u(I)' * 1/(2*area(j))*diag(signum)*diag(n)...   % 26/8/2
      *(ones(3,1)* coord(1,:) - coord);
   p2 = u(I)' * 1/(2*area(j))*diag(signum)*diag(n)...
      *(ones(3,1)* coord(2,:) -coord);
   p3 = u(I)' * 1/(2*area(j))*diag(signum)*diag(n)...
      *(ones(3,1)* coord(3,:) -coord);
   %p(m:m+2,:)=[p1;p2;p3];
   p(3*(j-1)+[1,2,3],:)=[p1;p2;p3]; %CB 9/4/3
   %m=m+3; 
   
   % =======================
   %% second way              % this gives the same result as the first way
   coord1=coordinate(element(j,:),:)';   
   N=coord1(:)*ones(1,3)-repmat(coord1,3,1);
   PC=reshape( N* 1/(2*area(j))*diag(signum)*diag(n)*u(I),2,3);
   %P(m:m+2,:)=[PC(1,:)',PC(2,:)'];
   P(3*(j-1)+[1,2,3],:)=[PC(1,:)',PC(2,:)'];%CB 9/4/3
   %m=m+3;     
end   
