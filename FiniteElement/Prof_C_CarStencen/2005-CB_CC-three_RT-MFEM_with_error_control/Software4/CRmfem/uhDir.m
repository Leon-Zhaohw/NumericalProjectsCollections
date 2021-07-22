
% C.Bahriawati
% Compute uh on T laying on \Gamma_D
% ==================================

function [uhD,calculated]=uhDir(coordinate,element,Dirichlet,nodes2edge,...
                          nodes2element,edge2element,ph,calculated,u)
                      
MidPoint=reshape(sum(reshape(coordinate(element',:),3,...
              2*size(element,1))),size(element,1),2)/3;

uhD=zeros(size(element,1),1);   
% ===========================================
% To compute u_h on T laying on Dir boundary

for k=1:size(Dirichlet,1)
   ElDir=nodes2element(Dirichlet(k,1),Dirichlet(k,2)) ;
   I=diag(nodes2edge(element(ElDir,[2 3 1]),element(ElDir,[3 1 2])));
   [ps,val]=DivPsi(coordinate,element,nodes2edge,edge2element,ElDir);
   %uI=([-1,1,1;1,-1,1;1,1,-1]*u(I));
   %valU=uI'*val;
   valU= ([-1,1,1;1,-1,1;1,1,-1]*u(I))'*val;
   uTimesPsi=[valU(1,1)*eye(3)  valU(1,2)*eye(3)  valU(1,3)*eye(3)]*ps;
   coordEval=[4 1 1;1 4 1;1 1 4]/6 * coordinate(element(ElDir,:),:);
   fbar=(det([1 1 1;coordinate(element(ElDir,:),:)'])* ...
             f(sum(coordEval )/2)) * ...
              (coordEval - ones(size( coordEval,1),1)*MidPoint(ElDir,:))/2;
   %fEval=fbar.*ps(1:3,:)+fbar.*ps(4:6,:)+fbar.*ps(7:9,:);
   %repmat(fbar,1,3)*ps
   %reshape(repmat(fbar,3,1).*ps,2,3)
   
   p=uTimesPsi - ( fbar.*ps(1:3,:)+fbar.*ps(4:6,:)+fbar.*ps(7:9,:) );
   uhD(ElDir)=det([1,1,1;coordinate(element(ElDir,:),:)']) * ...
              (sum( (  p)' ) * ones(3,1)/6 );
   calculated(ElDir) = 1;   
end


