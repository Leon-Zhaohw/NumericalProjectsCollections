global nu lambda mu H1 H sigma_y dt C1 C2 C3 C4 C5 C6 e0 sigma0 a0
th=1.0;time=[0:0.01:1];nu=0;
%initialization
load coordinates.dat;load elements.dat;load dirichlet.dat;load neumann.dat;
N=size(coordinates,1);initvector=zeros(2*N,1);NJ=size(elements,1);
sigma0=zeros(NJ,4);a0=zeros(NJ,1);u0=zeros(2*N,1);
%material parameters
lambda=4.142857142857144e+03;mu=1.035714285714286e+03;
C1=lambda+mu;C5=2*mu;sigma_y=0.1;H=1000;H1=1;
for step=2:length(time)
 t0=time(step-1);t1=time(step);dt=t1-t0;
 C2=1/(2*mu)*nu*(1+H^2*sigma_y^2)+th*dt*(1+H1*H^2*sigma_y^2/(2*mu));
 C4=H1*H^2*th*dt*sigma_y^2+nu*(1+H^2*sigma_y^2);
 C3=th*dt*sigma_y*(ones(NJ,1)+a0*H);C6=1/(2*mu)*(ones(NJ,1)+a0*H)*sigma_y; 
 e0=1/(4*(lambda+mu))*tr2(sigma0)*[1,0,0,1]+1/(2*mu)*dev2(sigma0);
 u_th=fem(coordinates,elements,dirichlet,neumann,initvector,u0,t0,t1,th,N);
 [sigma_th,a1_th]=tension(coordinates,elements,u_th,u0);
 sigma1=(sigma_th+(th-1)*sigma0)/th;a1=(a1_th+(th-1)*a0)/th;
 u1=(u_th+(th-1)*u0)/th;u0=u1;sigma0=sigma1;a0=a1;
 show(coordinates,elements,sigma1,u1,lambda,mu);
end










