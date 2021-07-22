function [M,F]=stima(M0,F0,loccoordinates,u1,u0,j);
global mu sigma_y C1 C2 C3 C4 e0
Dphi=inv([1,1,1,1;loccoordinates'])*[zeros(1,3);eye(3)];
eps=zeros(12,9);Deta1=zeros(12,9);Deta2=zeros(12,9);
Deta1(1:3:10,1:3)=Dphi;Deta1(2:3:11,4:6)=Dphi;Deta1(3:3:12,7:9)=Dphi;
Deta2(1:3:10,[1,4,7])=Dphi;Deta2(2:3:11,[2,5,8])=Dphi;Deta2(3:3:12,[3,6,9])=Dphi;
eps=(Deta1+Deta2)/2;v=(u1-u0)'*eps+e0(j,:);
T=det([1,1,1,1;loccoordinates'])/6;
if  norm(dev3(v))-sigma_y/(2*mu)>0
 C5=C2+C3/norm(dev3(v));C6=C3/norm(dev3(v))^3*(dev3(eps)*dev3(v)');
else
 C5=2*mu;C6=zeros(12,1);
end
M=M0+T*(C1*tr3(eps)*tr3(eps)'+C5*dev3(eps)*eps'-C6*dev3(v)*eps');
F=F0+T*(C1*tr3(v)*tr3(eps)+C5*eps*dev3(v)');
