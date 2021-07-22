% Function ph_OnRTElement.m computes $ph \in RT_0$
% from the non-conforming solution $uNC$ in the paper
% "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
% MFEM with a Posteriori Error Control" by C.Bahriawati and C. Carstensen


function ph=ph_OnRTElement(element,coordinate,nodes2edge,noedges,...
                           edge2element,uNC)
MidPoint=reshape(sum(reshape(coordinate(element',:),3,...
                             2*size(element,1))),size(element,1),2)/3;
ph=zeros(3*size(element,1),2);
for j=1:size(element,1)
  I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
  gradUNC=([-1,1,1;1,-1,1;1,1,-1]*uNC(I))'*...
          ([1,1,1;coordinate(element(j,:),:)']\[0,0;1,0;0,1]);
  ph(3*(j-1)+[1,2,3],:)=ones(3,1)*gradUNC-(det([1 1 1;...
  coordinate(element(j,:),:)'])*f(sum(coordinate(element(j,:),:))/2))*...
  (coordinate(element(j,:),:)-ones(3,1)*MidPoint(j,:))/2;
end