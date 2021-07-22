 global lambda mu sigma_y C1 C2 C3 C4 e0 
 th=1.0;time=[0:0.05:0.4];nu=0;
 %initialization
 load coordinates.dat; load elements.dat;load dirichlet.dat;load neumann.dat;
 N=size(coordinates,1);initvector=zeros(2*N,1);
 sigma0=zeros(size(elements,1),4);u0=zeros(2*N,1);
 %material parameters
 lambda=1.107438169066076e+05;mu=8.019379844961240e+04;C1=lambda+mu;sigma_y=450;
 for step=2:length(time)
  t0=time(step-1);t1=time(step);dt=t1-t0;
  C2=nu/(nu/(2*mu)+th*dt);C3=dt*sigma_y/(nu/(2*mu)+th*dt);
  e0=1/(4*(lambda+mu))*tr2(sigma0)*[1,0,0,1]+1/(2*mu)*dev2(sigma0);
  u_th=fem(coordinates,elements,dirichlet,neumann,initvector,u0,t0,t1,th,N);
  sigma_th=tension(coordinates,elements,u_th,u0,sigma0);
  u1=(u_th+(th-1)*u0)/th;sigma1=(sigma_th+(th-1)*sigma0)/th; 
  show(coordinates,elements,sigma1,u1,lambda,mu);
  u0=u1;sigma0=sigma1;
 end










