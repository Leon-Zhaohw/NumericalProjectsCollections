
% C.Bahriawati
% Search and compute uh on neighboring elements T
% ==============================================

function [uh,calculated]=uhN(element,coordinate,nodes2edge,edge2element,...
                          noedges,nodes2element,ph,uh,IntPh,calculated)

 weightedpoint=reshape(sum(reshape(coordinate(element',:),3,...
               2*size(element,1))),size(element,1),2)/3;                     
                      
listofelement=1:size(element,1);
neighbours=sparse(size(element,1),3);

for j=1:3
 neighbours(:,j) = diag(nodes2element(element(:,rem(j,3)+1),...
                                      element(:,rem(j-1,3)+1)));
end      

ell=size(element,1);
j = 1;
while j<=ell 
  if calculated(listofelement(j))
     j = j+1;
  elseif nnz(calculated(nonzeros(neighbours(listofelement(j),:))))
   % Compute uh on T+ and T-
   uh=computeUh(element,coordinate,nodes2edge,edge2element,...
                noedges,nodes2element,uh,IntPh,...
                listofelement(j),neighbours(listofelement(j),:));             
   calculated(listofelement(j))==1;
   [I,J]=find(neighbours==listofelement(j));
   for k=1:size(I,1)
     K=find(listofelement==I(k));
     if K>ell
     listofelement([ell+1,K])=listofelement([K,ell+1]);
     ell=ell+1;
     end
   end % of for
   j=j+1;
 else
   listofelement([j,ell])=listofelement([ell,j]);
   ell=ell-1;
 end  % of if
end
