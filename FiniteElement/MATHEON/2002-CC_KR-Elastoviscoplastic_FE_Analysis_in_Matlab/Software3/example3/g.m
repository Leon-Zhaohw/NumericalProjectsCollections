function sforce=g(x,n,t)
[phi,r]=cart2pol(x(:,1),x(:,2));
phiNeg =find(phi<0);
phi(phiNeg)=phi(phiNeg)+2*pi*ones(size(phiNeg));
for i= 1: size(x,1)
 if r(i)>.5 & r(i)<1.5
  sforce(i,:) = t*[cos(phi(i)),sin(phi(i))];
 elseif r(i)>1.5 & r(i)<2.5
  sforce(i,:)=-t/4*[cos(phi(i)),sin(phi(i))];
 end
 if x(i,1)==0 | x(i,2)==0
  sforce(i,:) = [0 0];
 end
end
