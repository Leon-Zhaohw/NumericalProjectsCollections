function [f,df,x]=rsgeno(N,rL,kl,ku,us,seed);
%rsgeno generates 1D random rough surfaces with bandlimited ocean spectrum.
%
%   [f,df,x]=rsgeno(N,rL,kl,ku,us,seed)
%
%   INPUT:
%
%   N=total number of sample points
%   rL=rough surface length
%   kl=lower wavenumber cutoff
%   ku=upper wavenumber cutoff
%   us=wind friction velocity
%	 seed=seed of random number generator
%
%   OUTPUT:
%
%   f=rough surface profile
%   df=df/dx
%   x=sample points on the surface
%
% -- Part of the Electromagnetic Wave MATLAB Library (EWML) --
%    <http://www.emwave.com/>

% Original: C. O. Ao, 2000.

randn('seed',seed);
y=randn(N,1);
for n=1:(N/2-1);
  bh(n)=(y(2*n-1)+i*y(2*n))/sqrt(2);
end;

bhc=conj(bh);
bhf=fliplr(bhc);
bi=[bh y(N-1) bhf y(N)];
kx=2.*pi*[-N/2+1:1:N/2]/rL;

y1=zeros(size(kx));
i1=find(abs(kx) >= kl & abs(kx) <= ku);
y1(i1)=sqrt(wkos(kx(i1),us));
y=y1*sqrt(2*pi*rL);
b=y.*bi;
xs=[b(N/2+1:1:N) b(1:1:N/2)];
xt=[xs(N),xs(1:1:N-1)];
ft=ifft(xt,N);
ft=ft*N/rL;
fs=[ft(2:1:N),ft(1)];
f=[fs(N/2+1:1:N) fs(1:1:N/2)];
f=real(f);
dx=rL/N;
x=[-N/2+1:1:N/2]*dx;
n=2:N-1;
df1=(f(n+1)-f(n-1))/(2*dx);
df=[(f(2)-f(N))/(2*dx),df1,(f(1)-f(N-1))/(2*dx)];

%%%%%%%%%%%%%%%%%%%%%
% 1D ocean spectrum %
%%%%%%%%%%%%%%%%%%%%%
function y=wkos(kx,us)

kr=abs(kx);
%a0=0.008;
a0=0.004;
a=0.225;
b=1.25;
kj=2;
gs=9.81+7.25e-5*kx.^2;
z0=6.84e-5/us+4.28e-3*us^2-4.43e-4;
U195=us/0.4*log(19.5/z0);
kc=9.81/U195^2;

i1=find(kr>kj);
kr1=kr(i1);
y(i1)=(a0./kr1.^3).*(b*kr1*us^2./gs(i1)).^(a*log10(kr1/kj));

i2=find(kr<=kj);
kr2=kr(i2);
y(i2)=(a0./kr2.^3).*exp(-0.74*(kc./kr2).^2);
