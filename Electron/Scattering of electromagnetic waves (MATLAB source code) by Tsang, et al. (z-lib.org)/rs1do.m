function [tsd,sig]=rs1do(wave,nr,N,rL,kl,ku,us,g,tid,nsa,seed)
%rs1do computes the bistatic scattering coefficient for Gaussian rough surfaces 
%   with ocean spectrum.
%
%   [tsd,sig]=rs1do(wave,nr,N,rL,kl,ku,us,g,tid,nsa,seed)
%
%   INPUT:
%
%   wave=wavelength
%   nr=total number of surface realizations
%   N=total number of sample points
%   rL=rough surface length
%   kl=lower wavenumber cutoff
%   ku=upper wavenumber cutoff
%   us=wind friction velocity
%   g=tapering parameter for incident wave
%   tid=incident angle in degree
%   nsa=number of scattered angles from -90 deg to 90 deg
%   seed=seed for random number generator
%
%   OUTPUT:
%
%   tsd=scattered angles in degree
%   sig=bistatic scattering coef (MoM)
%
%	 REQUIRES: rsgeno.m for generation of rough surfaces with Ocean Spectrum
%
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: C. O. Ao, 2000

tai=tid*pi/180;
k=2*pi/wave;
ti=tan(tai);
ci=cos(tai);
denom=1+2*ti^2;
denom=denom/(2*(k*g*ci)^2);
denom=8*pi*k*g*sqrt(pi/2)*ci*(1-denom);

euler=1.78107;
e=exp(1);

sig=zeros(nsa,1);
randn('seed',seed);
for ir=1:nr
  fprintf('Processing realization %i ...\n',ir);
  [f,df,x]=rsgeno(N,rL,kl,ku,us,seed);
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