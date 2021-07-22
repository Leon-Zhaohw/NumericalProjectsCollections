% C.Bahriawati
% Compute Integral of ph.psi_E on T 
% ==================================


function IntPh=IntPhOmega(element,coordinate,u,nodes2edge,edge2element)

MidPoint = reshape(sum(reshape(coordinate(element',:),3,...
                2*size(element,1))),size(element,1),2)/3;
 
IntPh=0;
for j = 1:size(element,1)
  % area(j) = det([1 1 1;coordinate(element(j,:),:)'])/2;
  % coord=coordinate(element(j,:),:) ;
   I=diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
   [ps,val]=DivPsi(coordinate,element,nodes2edge,edge2element,j);
   valU= ([-1,1,1;1,-1,1;1,1,-1]*u(I))'*val;
   uTimesPsi=[valU(1,1)*eye(3)  valU(1,2)*eye(3)  valU(1,3)*eye(3)]*ps;
   coordEval=[4 1 1;1 4 1;1 1 4]/6 * coordinate(element(j,:),:);
   fbar=(det([1 1 1;coordinate(element(j,:),:)'])* ...
             f(sum(coordEval )/2)) * ...
              (coordEval - ones(size(coordEval,1),1)*MidPoint(j,:))/2;
          
   %fEval=fbar.*ps(1:3,:)+fbar.*ps(4:6,:)+fbar.*ps(7:9,:);
   PhPsi= uTimesPsi - ( fbar.*ps(1:3,:)+fbar.*ps(4:6,:)+fbar.*ps(7:9,:));
   IntPh= IntPh + det([1,1,1;coordinate(element(j,:),:)']) * ...
                  (sum((PhPsi)') * ones(3,1)/6);
end

