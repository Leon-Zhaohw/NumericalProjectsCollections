function [tsd,sig,sigka,sigspm]=rs1dg(wave,nr,N,rL,h,lc,g,tid,nsa,seed)
%rs1dg computes the bistatic scattering coefficient for Gaussian rough surfaces 
%   with Gaussian spectrum.
%
%   [tsd,sig,sigka,sigspm]=rs1dg(wave,nr,N,rL,h,lc,g,tid,nsa,seed)
%
%   INPUT:
%
%   wave=wavelength
%   nr=total number of surface realizations
%   N=total number of sample points
%   rL=rough surface length
%   h=rms height
%   lc=correlation length
%   g=tapering parameter for incident wave
%   tid=incident angle in degree
%   nsa=number of scattered angles from -90 deg to 90 deg
%   seed=seed for random number generator
%
%   OUTPUT:
%
%   tsd=scattered angles in degree
%   sig=bistatic scattering coef (MoM)
%   sigka=bistatic scattering coef (KA)
%   sigspm=bistatic scattering coef (SPM)
%
%	 REQUIRES: rsgeng.m for generation of Gaussian rough surfaces
%
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: L. Tsang, 1998

tai=tid*pi/180;
h2=h^2;
k=2*pi/wave;
ti=tan(tai);
ci=cos(tai);
denom=1+2*ti^2;
denom=denom/(2*(k*g*ci)^2);
denom=8*pi*k*g*sqrt(pi/2)*ci*(1-denom);

sig=zeros(nsa,1);
randn('seed',seed);
for ir=1:nr
  fprintf('Processing realization %i ...\n',ir);
  [f,df,x]=rsgeng(N,rL,h,lc,seed);
  seed=randn('seed');
  b=incid(k,x,f,tai,g);
  dx=rL/N;
  [xm,xn]=meshgrid(x);
  [fm,fn]=meshgrid(f);
  arg=k*sqrt((xm-xn).^2+(fm-fn).^2);
  for m=1:N
	  arg(m,m)=1;
  end
  A=dx*i*besselh(0,1,arg)/4;
  euler=1.78107;
  e=exp(1);
  for m=1:N
    dl=sqrt(1+df(m)*df(m))*dx;
    A(m,m)=(i/4)*dx*(1+(2*i/pi)*(log(euler*k*dl/4)-1));
  end
  b=b.';
  u=A\b;
  u=u.';
  
  dan=180/(nsa+1);
  for m=1:nsa
    tsd(m)=-90+m*dan;
    tas=tsd(m)*pi/180;
    ss=sin(tas);
    cs=cos(tas);
    integ=exp(-i*k*(ss*x+f*cs));
    integ=integ.*u*dx;
    psis=sum(integ);
    sige=abs(psis)^2/denom;
    sig(m)=(sige+(ir-1)*sig(m))/ir;
  end
end

% Kirchhoff Approximation
for m=1:nsa;
  tas=tsd(m)*pi/180;
  ss=sin(tas);
  cs=cos(tas);
  csa(m)=cs^2;
  kxq=k*(ss-sin(tai));
  kx(m)=kxq;
  csum=ci+cs;
  fac=(1+cos(tai+tas))^2*k^3/ci*exp(-(k*h*csum)^2);
  p=0;
  arr=1/(k*csum)^2;
  for ma=1:20;
    term=lc/(2*sqrt(pi));
    term=term/sqrt(ma);
    term=term*exp(-(kxq*lc/2)^2/ma);
    arr=arr*(k*csum*h)^2/ma;
    term1=arr*term;
    p=p+term1;
  end
  sigka(m)=p*fac;
end;

%Small Perturbation Method
sigspm=csa.*wk(kx,h,lc)*4*k^3*cos(tai);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              incid.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b=incid(k,x,z,tai,g)
% generates the spatial tapered incident wave

ti=tan(tai);
ci=cos(tai);
si=sin(tai);
fac=((x+z*ti)/g).^2;
kg=(k*g*ci)^2;
w=(2*fac-1)/kg;
b=exp(i*k*(x*si-z*ci).*(1+w)-fac);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               wk.m                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=wk(kx,h,L)
% Gaussian spectral density

y=h^2*L*exp(-(kx*L*0.5).^2)/(2*sqrt(pi));
