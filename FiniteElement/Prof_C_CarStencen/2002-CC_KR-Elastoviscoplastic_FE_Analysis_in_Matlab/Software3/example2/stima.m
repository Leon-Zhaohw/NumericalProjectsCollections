function [M,F]=stima(M0,F0,lcoordinates,u1,u0,j);
global C1 C2 C3 C4 C5 C6 e0
Dphi=inv([1,1,1;lcoordinates'])*[0,0;eye(2)];
eps=zeros(6,4);Deta1=zeros(6,4);Deta2=zeros(6,4);
Deta1(1:2:5,1:2)=Dphi;Deta1(2:2:6,3:4)=Dphi;
Deta2(1:2:5,[1,3])=Dphi;Deta2(2:2:6,[2,4])=Dphi;
eps=(Deta1+Deta2)/2;v=(u1-u0)'*eps+e0(j,:);
T=det([1,1,1;lcoordinates'])/2;
if norm(dev2(v))-C6(j)>0
 C7=C3(j)/(C2*norm(dev2(v)))+C4/C2;C8=C3(j)/(C2*norm(dev2(v))^3)*dev2(eps)*dev2(v)';
else    
  C7=C5;C8=zeros(6,1);
end
M=M0+T*(C1*tr2(eps)*tr2(eps)'+C7*dev2(eps)*eps'+C8*dev2(v)*eps');
F=F0+T*(C1*tr2(v)*tr2(eps)+C7*eps*dev2(v)');
  
 
 
     
 














