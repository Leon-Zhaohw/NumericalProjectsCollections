function [sigma1,a1]=tension(coordinates,elements,u_1,u_0);
global H1 H nu sigma_y dt C1 C2 C3 C4 C5 C6 e0 a0 sigma0
sigma1=sigma0; a1=a0;
for j=1:size(elements,1)
 lcoordinates=coordinates(elements(j,:),:);
 Dphi=inv([1,1,1;lcoordinates'])*[0,0;eye(2)];
 Deta1=zeros(6,4);Deta2=zeros(6,4);u1=zeros(6,1);u0=zeros(6,1);
 Deta1(1:2:5,1:2)=Dphi;Deta1(2:2:6,3:4)=Dphi;
 Deta2(1:2:5,[1,3])=Dphi;Deta2(2:2:6,[2,4])=Dphi;
 u1(1:2:5)=u_1(2*elements(j,:)-1);u1(2:2:6)=u_1(2*elements(j,:));
 u0(1:2:5)=u_0(2*elements(j,:)-1);u0(2:2:6)=u_0(2*elements(j,:));
 eps=(Deta1+Deta2)/2;v=(u1-u0)'*eps+e0(j,:);
 if norm(dev2(v))-C6(j)>0
  sigma1(j,:)=C1*tr2(v)*[1,0,0,1]+(C3(j)/(C2*norm(dev2(v)))+C4/C2)*dev2(v);
  a1(j)=(H1*H*dt*sigma_y*(norm(dev2(sigma1(j,:)))-sigma_y)+...
         a0(j)*nu*(1+H^2*sigma_y^2))/(nu*(1+H^2*sigma_y^2)+dt*H1*H^2*sigma_y^2);
 else
  sigma1(j,:)=C1*tr2(v)*[1,0,0,1]+C5*dev2(v);
 end
end







