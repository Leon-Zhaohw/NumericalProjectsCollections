function [M,F]=STEMA(M0,F0,loccoordinates,u1,u0,j);
global mu sigma_y C1 C2 C3 C4 e0 al0
Dphi=inv([1,1,1;loccoordinates'])*[0,0;eye(2)];
eps=zeros(6,4);Deta1=zeros(6,4);Deta2=zeros(6,4);
Deta1(1:2:5,1:2)=Dphi;Deta1(2:2:6,3:4)=Dphi;
Deta2(1:2:5,[1,3])=Dphi;Deta2(2:2:6,[2,4])=Dphi;
eps=(Deta1+Deta2)/2;v0=(u1-u0)'*eps+e0(j,:);v1=-1/(2*mu)*al0(j,:);v=v0+v1;
T=det([1,1,1;loccoordinates'])/2;
if norm(dev2(v))-sigma_y/(2*mu)>0
 C5=C2+C3/norm(dev2(v));C6=(C3/norm(dev2(v))^3)*dev2(eps)*dev2(v)';d=1;   
else
 C5=C4;C6=zeros(6,1);d=0;
end
M=M0+T*(C1*tr2(eps)*tr2(eps)'+C5*dev2(eps)*eps'-C6*dev2(v)*eps');
F=F0+T*(C1*tr2(v0)*tr2(eps)+C5*eps*dev2(v+d*v1)'+d*eps*dev2(al0(j,:))');
  
     
     
     
     
 















