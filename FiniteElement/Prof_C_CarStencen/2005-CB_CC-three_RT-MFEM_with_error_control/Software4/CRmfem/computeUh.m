
 %nEl=listofelement(j);  % CB
 %NEl=neighbours(listofelement(j),:);                     

function uh=computeUh(element,coordinate,nodes2edge,edge2element,...
                      noedges,nodes2element,uh,IntPh,nEl,NEl); 
                      
 
[ps,val]=DivPsi(coordinate,element,nodes2edge,edge2element,nEl);
valInt=val(1,1)*ps(1:3,:) + val(2,2)*ps(4:6,:) + val(3,3)*ps(7:9,:);
val1=uh(nEl)* det([1,1,1;coordinate(element(nEl,:),:)']) * ...
     (sum( ( valInt)' ) * ones(3,1)/6 );
for k=1:size(NEl,1)
   N=NEl(k);
   [ps,val]=DivPsi(coordinate,element,nodes2edge,edge2element,N);
   valInt=[val(1,1)*eye(3)  val(1,2)*eye(3)  val(1,3)*eye(3)]*ps;
   IntDiv=det([1,1,1;coordinate(element(NEl(k),:),:)']) * ...
          (sum( (valInt)' ) * ones(3,1)/6 );
   uh(NEl(k))= (-IntPh-val1)/IntDiv;
end 