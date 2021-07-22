function ut=FEM3_d(coordinates,elements,dirichlet,neumann,ut,u0,t0,t1,th,N,NJ);
%Initialisation
DF=sparse(3*N,3*N);
Q=zeros(3*N,1);P=zeros(3*N,1);
list=3*elements(:,[1,1,1,2,2,2,3,3,3,4,4,4])-ones(NJ,1)*[2,1,0,2,1,0,2,1,0,2,1,0]; 
%Assembly
for j=1:size(elements,1)
 L=list(j,:);
 [DF(L,L),Q(L)]=stima(DF(L,L),Q(L),coordinates(elements(j,:),:),ut(L),u0(L),j);
end
%Volume forces
for j=1:size(elements,1)
 lcoordinates=coordinates(elements(j,:),:);vol=det([1,1,1,1;lcoordinates'])/6;
  b=vol*(1-th)*f(sum(lcoordinates,1)/3,t0)'/4+th*f(sum(lcoordinates,1)/3,t1)'/4;
  P(list(j,:))=P(list(j,:))+repmat(b,4,1); 
 end
%Neumann conditions
 if ~isempty(neumann)
  Nlist=3*neumann(:,[1,1,1,2,2,2,3,3,3])-ones(size(neumann,1),1)*repmat([2,1,0],1,3);
  for j=1:size(neumann,1)
   NV=cross(coordinates(neumann(j,2),:)-coordinates(neumann(j,1),:), ...
            coordinates(neumann(j,3),:)-coordinates(neumann(j,1),:));
   area=norm(NV)/2;NV=NV/norm(NV);
   gs=area*((1-th)*g(sum(coordinates(neumann(j,:),:))/2,NV,t0)'/2 ...
      +th*g(sum(coordinates(neumann(j,:),:))/2,NV,t1)'/2);
   P(Nlist(j,:))=P(Nlist(j,:))+repmat(gs,3,1);
  end
 end 
 F=Q-P;
%Dirichlet conditions
 dnodes=unique(dirichlet);
 [W0,M0]=u_D(coordinates(dnodes,:),t0);[W1,M1]=u_D(coordinates(dnodes,:),t1);
 W=(1-th)*W0+th*W1;M=(1-th)*M0+th*M1;B=sparse(size(W,1),3*N);
 for k = 0:2
  for j = 0:2
   B(1+j:3:size(M,1),3*dnodes-2+k)=diag(M(1+j:3:size(M,1),1+k));
  end
 end
 mask=find(sum(abs(B)'));freeNodes=find(~sum(abs(B)));  
 B=B(mask,:);W=W(mask,:);norm0=norm(F(freeNodes));nit=0;
 while (norm(F(freeNodes))>10^(-10)+10^(-6)*norm0) & (nit<100)
  nit=nit+1;
% Calculating the solution
  F=[DF*ut-F;W];DF=[DF,B';B,sparse(size(B,1),size(B,1))];  
  x=DF\F;ut=x(1:3*N);lb=x(3*N+1:size(x,1));
% Assembly 
  DF=sparse(3*N,3*N);Q=zeros(3*N,1);
  for j=1:size(elements,1)
   L=list(j,:);
   [DF(L,L),Q(L)]=stima(DF(L,L),Q(L),coordinates(elements(j,:),:),ut(L),u0(L),j);
  end
  F=Q-P;
 end
 if (nit==100)  disp('Fails to converge within 100 iteration steps!'); end  

 
 
  






