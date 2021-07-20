function sigma1=Spannung1(coordinates,elements,u_1,u_0,sigma0);
global mu sigma_y C1 C2 C3 e0
sigma1=sigma0;
for j=1:size(elements,1)
 lcoordinates=coordinates(elements(j,:),:);
 Dphi=inv([1,1,1,1;lcoordinates'])*[zeros(1,3);eye(3)];
 eps=zeros(12,9);Deta1=zeros(12,9);Deta2=zeros(12,9);
 Deta1(1:3:10,1:3)=Dphi;Deta1(2:3:11,4:6)=Dphi;Deta1(3:3:12,7:9)=Dphi;
 Deta2(1:3:10,[1,4,7])=Dphi;Deta2(2:3:11,[2,5,8])=Dphi;Deta2(3:3:12,[3,6,9])=Dphi;
 u1=u_1(3*elements(j,reshape(repmat([1:4],3,1),1,12))-repmat([2,1,0],1,4));
 u0=u_0(3*elements(j,reshape(repmat([1:4],3,1),1,12))-repmat([2,1,0],1,4));
 eps=(Deta1+Deta2)/2;v=(u1-u0)'*eps+e0(j,:);
 if norm(dev3(v))-sigma_y/(2*mu)>0
  sigma1(j,:)=C1*tr3(v)*[1,0,0,0,1,0,0,0,1]+(C2+C3/norm(dev3(v)))*dev3(v);   
 else
  sigma1(j,:)=C1*tr3(v)*[1,0,0,0,1,0,0,0,1]+2*mu*dev3(v);
 end
end








