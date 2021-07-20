%function tmp=signedge(edge2element,x)
function tmp=signedge(r)
   % check sign
   
   if ~isempty(r)  
     tmp=zeros(3,1);
     tmp(r,1)=-1;
     k=find(~tmp);
     tmp(k,1)=1; 
   else
     tmp=ones(3,1);
   end