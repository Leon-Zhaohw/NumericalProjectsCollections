 global lambda mu nu sigma_y dt C1 C2 C3 C4 k1 e0 sigma0 al0
 th=1.0;time=[0:0.01:1];nu=0;
%initialization
 load coordinates.dat;load elements.dat;load dirichlet.dat;load neumann.dat;
 N=size(coordinates,1);initvector=zeros(2*N,1);NJ=size(elements,1);   
 sigma0=zeros(NJ,4);al0=zeros(NJ,4);u0=zeros(2*N,1);
 %material parameters
 lambda=5.108359133126935e+04;mu=2.631578947368421e+04;
 C1=lambda+mu;C4=2*mu;sigma_y=10^(-8);k1=1;
 for step=2:length(time)
  t0=time(step-1);t1=time(step);dt=t1-t0;
  C2=(th*dt*k1+2*nu)/(th*dt+1/(2*mu)*th*dt*k1+nu/mu);
  C3=th*dt*sigma_y/(th*dt+1/(2*mu)*th*dt*k1+nu/mu);
  e0=1/(4*(lambda+mu))*tr2(sigma0)*[1,0,0,1]+1/(2*mu)*dev2(sigma0);
  u_th=fem(coordinates,elements,dirichlet,neumann,initvector,u0,t0,t1,th,N);
  [sigma_th,al1_th]=tension(coordinates,elements,u_th,u0);
  sigma1=(sigma_th+(th-1)*sigma0)/th ; 
  al1=(al1_th+(th-1)*al0)/th;
  u1=(u_th+(th-1)*u0)/th;u0=u1;sigma0=sigma1;al0=al1;
  show(coordinates,elements,sigma1,u1,lambda,mu);
end









