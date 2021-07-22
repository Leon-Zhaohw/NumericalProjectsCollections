function [sigma1,al1]=tension(coordinates,elements,u_1,u_0);
global lambda mu sigma_y C1 C2 C3 C4 k1 e0 sigma0 al0
sigma1=sigma0;al1=al0;
for j=1:size(elements,1)
 lcoordinates=coordinates(elements(j,:),:);
 Dphi=inv([1,1,1;lcoordinates'])*[0,0;eye(2)];
 Deta1=zeros(6,4);Deta2=zeros(6,4);u1=zeros(6,1);u0=zeros(6,1);
 Deta1(1:2:5,1:2)=Dphi;Deta1(2:2:6,3:4)=Dphi;
 Deta2(1:2:5,[1,3])=Dphi;Deta2(2:2:6,[2,4])=Dphi;
 u1(1:2:5)=u_1(2*elements(j,:)-1);u1(2:2:6)=u_1(2*elements(j,:));
 u0(1:2:5)=u_0(2*elements(j,:)-1);u0(2:2:6)=u_0(2*elements(j,:));
 eps=(Deta1+Deta2)/2;v=(u1-u0)'*eps+e0(j,:);
 if norm(dev2(v-1/(2*mu)*al0(j,:)))-1/(2*mu)*sigma_y>0
  sigma1(j,:)=C1*tr2(v)*[1,0,0,1]+(C2+C3/norm(dev2(v-1/(2*mu)*al0(j,:))))*...
              dev2(v-1/(2*mu)*al0(j,:))+dev2(al0(j,:));
  al1(j,:)=al0(j,:)+k1*(v-1/(4*(lambda+mu))*tr2(sigma1(j,:))*[1,0,0,1]-...
           1/(2*mu)*dev2(sigma1(j,:)));
 else
  sigma1(j,:)=C1*tr2(v)*[1,0,0,1]+C4*dev2(v);al1(j,:)=al0(j,:);
 end
end

