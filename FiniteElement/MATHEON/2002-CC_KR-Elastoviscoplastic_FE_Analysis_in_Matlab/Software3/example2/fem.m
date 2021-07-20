function ut=FEM(coordinates,elements,dirichlet,neumann,ut,u0,t0,t1,th,N);
%Initialisation
 DF=sparse(2*N,2*N);Q=zeros(2*N,1);P=zeros(2*N,1);
 list=2*elements(:,[1,1,2,2,3,3])-ones(size(elements,1),1)*[1,0,1,0,1,0];   
%Assembly
 for j=1:size(elements,1)
  L=list(j,:);
  [DF(L,L),Q(L)]=stima(DF(L,L),Q(L),coordinates(elements(j,:),:),ut(L),u0(L),j);
 end
%Volume forces
 for j=1:size(elements,1)
  lcoordinates=coordinates(elements(j,:),:);area=det([1,1,1;lcoordinates'])/2;
  b=area*(1-th)*f(sum(lcoordinates,1)/3,t0)'/3+th*f(sum(lcoordinates,1)/3,t1)'/3;
  P(list(j,:))=P(list(j,:))+repmat(b,3,1); 
 end
%Neumann conditions
 if ~isempty(neumann)
  Nlist=2*neumann(:,[1,1,2,2])-ones(size(neumann,1),1)*repmat([1,0],1,2);
  EV=coordinates(neumann(:,2),:)-coordinates(neumann(:,1),:);
  Lg=sqrt(sum(EV.*EV,2));EV=EV./[Lg,Lg];NV=EV*[0,-1;1,0];
  for j=1:size(neumann,1)
   gs=Lg(j)*((1-th)*g(sum(coordinates(neumann(j,:),:))/2,NV(j,:),t0)'/2 ...
       +th*g(sum(coordinates(neumann(j,:),:))/2,NV(j,:),t1)'/2);
   P(Nlist(j,:))=P(Nlist(j,:))+repmat(gs,2,1);
  end
 end 
 F=Q-P;
%Dirichlet conditions
 dnodes=unique(dirichlet);
 [W0,M0]=u_D(coordinates(dnodes,:),t0);[W1,M1]=u_D(coordinates(dnodes,:),t1);
 W=(1-th)*W0+th*W1;M=(1-th)*M0+th*M1;B=sparse(size(W,1),2*N);
 for k=0:1
  for j=0:1
   B(1+j:2:size(M,1),2*dnodes-1+k)=diag(M(1+j:2:size(M,1),1+k));
  end
 end
 maske=find(sum(abs(B)'));freeNodes=find(~sum(abs(B)));
 B=B(maske,:);W=W(maske,:);norm0=norm(F(freeNodes));nit=0;
 while (norm(F(freeNodes))>10^(-10)+10^(-6)*norm0) & (nit<100)
  nit=nit+1;
% calculation of the solution
  F=[DF*ut-F;W];DF=[DF,B';B,sparse(size(B,1),size(B,1))];  
  x=DF\F;ut=x(1:2*N);lb=x(2*N+1:size(x,1));
% Assembly
  DF=sparse(2*N,2*N);Q=zeros(2*N,1);
  for j=1:size(elements,1)
   L=list(j,:);
   [DF(L,L),Q(L)]=stima(DF(L,L),Q(L),coordinates(elements(j,:),:),ut(L),u0(L),j);
  end
  F=Q-P;
 end
 if (nit==100)  disp('Fails to converge within 100 iteration steps!'); end 

 
 
  






