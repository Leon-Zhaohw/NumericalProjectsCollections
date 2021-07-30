%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            ptscext.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main program for the Monte Carlo simulations of volumetric 
% scattering by clusters of point scatterers
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: by K.H. Ding, November 1998

% input parameters
nc=input('Enter number of clusters : ');
npc=input('Enter number of point scatterers per cluster : ');
wavel=input('Enter incident wavelenght (meter) : ');
L=input('Enter cubic box length (# of wavelengths) : ');
lc=input('Enter cluster length (# of wavelengths) : ');
fp=input('Enter real part of scattering amplitude : ');
tai=input('Enter polar incidence angle (degree) : ');
phi=input('Enter azimuthal incidence angle (degree) : ');
nrlz=input('Enter number of realizations : ');
seed=input('Enter seed for random numbers : ');

ntas=10;
gmatas=(1:ntas-1)./sqrt((2*(1:ntas-1)).^2-1);
[evtas,xatas]=eig(diag(gmatas,-1)+diag(gmatas,1));
[xatas,katas]=sort(diag(xatas));
watas=2*evtas(1,katas)'.^2;

nphs=20;
gmaphs=(1:nphs-1)./sqrt((2*(1:nphs-1)).^2-1);
[evphs,xaphs]=eig(diag(gmaphs,-1)+diag(gmaphs,1));
[xaphs,kaphs]=sort(diag(xaphs));
waphs=2*evphs(1,kaphs)'.^2;

dtr=pi/180;
twopi=2*pi;

N=nc*npc;
ncps=L/lc;
k=twopi/wavel;
fpp=twopi*fp*fp;
f=(fp+i*fpp)*wavel;
hsl=L*wavel/2;
csl=lc*wavel;
hcsl=csl/2;
V=(2*hsl)^3;
n0=N/V;
Kei=n0*4*pi*imag(f)/k;

% scattering directions
nsang=ntas*nphs;
kxs=ones(nsang,1);
kys=ones(nsang,1);
kzs=ones(nsang,1);
wf=ones(nsang,1);
iang=0;
for ita=1:ntas
  wt=watas(ita);
  tas=acos(xatas(ita));
  for iph=1:nphs
    iang=iang+1;
    wp=waphs(iph);
    wf(iang)=wt*wp;
    phs=pi+pi*xaphs(iph);
    cps=cos(phs);
    sps=sin(phs);
    cts=cos(tas);
    sts=sin(tas);
    kxs(iang)=-i*k*sts*cps;
    kys(iang)=-i*k*sts*sps;
    kzs(iang)=-i*k*cts;
    tdeg(iang)=tas/dtr;
    pdeg(iang)=phs/dtr;
  end
end

% incident direction
tai=tai*dtr;
phi=phi*dtr;
cti=cos(tai);
sti=sin(tai);
cpi=cos(phi);
spi=sin(phi);
kxi=i*k*sti*cpi;
kyi=i*k*sti*spi;
kzi=-i*k*cti;

fext=fopen('ptscext.dat','w+');

rand('seed',seed);

xrow=zeros(1,N);
yrow=zeros(1,N);
zrow=zeros(1,N);
xc=zeros(npc,nc);
yc=zeros(npc,nc);
zc=zeros(npc,nc);
xccol=zeros(N,1);
yccol=zeros(N,1);
zccol=zeros(N,1);
rx=zeros(N,N);
ry=zeros(N,N);
rz=zeros(N,N);
rr=ones(N,N);
rr1=ones(N,N);
b=ones(N,1);
Z=ones(N,N);
Psi=ones(N,1);
erow=ones(1,N);
ecol=ones(N,1);
emat=ones(npc,nc);
xcmat=eye(nc);
ycmat=eye(nc);
zcmat=eye(nc);

FF=ones(nsang,1);
FFc=zeros(nsang,1);
Pt=zeros(nsang,1);
Pc=zeros(nsang,1);
P=zeros(nsang,1);

% Monte Carlo Simulation
for ir=1:nrlz
  ir1=ir-1;

% positions of point scatterers
  if nc == 1
    xrow=hsl*(erow-2.0*rand(1,N));
    yrow=hsl*(erow-2.0*rand(1,N));
    zrow=hsl*(erow-2.0*rand(1,N));
  else
    xcmat=diag(hcsl+csl*floor(ncps*rand(1,nc)));
    ycmat=diag(hcsl+csl*floor(ncps*rand(1,nc)));
    zcmat=diag(hcsl+csl*floor(ncps*rand(1,nc)));
    xc=emat*xcmat;
    yc=emat*ycmat;
    zc=emat*zcmat;
    xccol=xc(:);
    yccol=yc(:);
    zccol=zc(:);
    xrow=xccol'+hcsl*(erow-2.0*rand(1,N));
    yrow=yccol'+hcsl*(erow-2.0*rand(1,N));
    zrow=zccol'+hcsl*(erow-2.0*rand(1,N));
  end

% incident
  b=exp(kxi*xrow'+kyi*yrow'+kzi*zrow');

% calculate separation between pairs of point scatterers
  rx=xrow'*erow-ecol*xrow;
  ry=yrow'*erow-ecol*yrow;
  rz=zrow'*erow-ecol*zrow;
  rr=sqrt(rx.^2+ry.^2+rz.^2);
  rr1=rr;
  for j=1:N
    rr1(j,j)=1.0;
  end

% impedance matrix
  Z=-f*exp(i*k*rr)./rr1;
  for j=1:N
    Z(j,j)=1.0;
  end

% exciting field
  Psi=Z\b;

% total scattering amplitude
  FF=f*exp(kxs*xrow+kys*yrow+kzs*zrow)*Psi;
  Pt=(ir1*Pt+abs(FF).*abs(FF))/ir;

% cohenent scattering amplitude
  FFc=(ir1*FFc+FF)/ir;
  Pc=abs(FFc).*abs(FFc);

% incoherent scattered power
  P=Pt-Pc;

% extinction
  Kec=pi*wf'*P/V;
  fprintf(fext,'%5u %9.6e %9.6e %9.6e \n',ir,Kec*wavel,Kei*wavel,Kec/Kei);
  fprintf('\n     realization = %6u \n',ir);
end
% N realizations done
fclose(fext);

% plot
clear;
load ptscext.dat -ascii;
figure;
plot(ptscext(:,1),ptscext(:,4));
xlabel('number of realizations ');
ylabel('ext. coeff. (cluster)/(uniform)');
title('      Convergence Test');
clear;