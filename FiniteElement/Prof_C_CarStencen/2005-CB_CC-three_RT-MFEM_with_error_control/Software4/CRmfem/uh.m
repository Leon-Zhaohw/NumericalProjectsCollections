
function uh=uhN(element,coordinate,nodes2edge,edge2element,noedges,NEl,nEl,pEval); 

weightedpoint = reshape(sum(reshape(coordinate(element',:),3,...
                2*size(element,1))),size(element,1),2)/3;
            
uh=zeros(size(element,1),1);
for j=1:3
   coord=coordinate(element(NEl(j),:),:);
   Area = det([1 1 1;coordinate(element(NEl(j),:),:)'])/2; 
   I=diag(nodes2edge(element(NEl(j),[2 3 1]),element(NEl(j),[3 1 2])));
   [r,s]=find(NEl(j)==edge2element(I,4));
   tmp=signedge(r);
   c=coordinate(element(NEl(j) ,[2 3 1]),:)-coordinate(element(NEl(j),[3 1 2]),:);
   n=[norm(c(1,[1 2])),norm(c(2,[1 2])),norm(c(3,[1 2]))];
   
  
   % ===================       
   % compute  \int_{\Omega}ph \psi_E dx
   psi = 1/(2*Area)*diag(tmp)*diag(n)...
      *(ones(3,1)*weightedpoint(NEl(j),:)-coord)
   p= pEval( NEl(j):NEl(j)+2,:);
   PhPsi = p'*psi;
   IntP = det([1,1,1;coordinate(element(ElDir,:),:)']) * ...
                         (PhPsi' * [4 1 1;1 4 1;1 1 4]/36)' ;  
   uh(NEL(j))= IntP-uh(nEL);
end
