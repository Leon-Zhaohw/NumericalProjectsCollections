function sforce=g(x,n,t)
sforce=zeros(size(x,1),2);
if (n(1)==1)
 sforce(:,2)=exp(2*t);
end
