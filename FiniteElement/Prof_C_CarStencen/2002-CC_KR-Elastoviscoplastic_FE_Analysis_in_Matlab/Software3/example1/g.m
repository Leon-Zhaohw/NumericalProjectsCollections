function sforce=g(x,n,t)
% surface force
sforce=zeros(size(x,1),2);
if (n(2)==1)
  sforce(:,2)=600*t ;
end
