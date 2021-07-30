function [f,df,x]=rsgenf(N,rL,h,kl,ku,Nf,sdim,seed);
%rsgenf generates bandlimited random rough surfaces using the 
%   Wierstrass-Mandelbrot function.
%
%   [f,df,x]=rsgenf(N,rL,h,kl,ku,Nf,sdim,seed)
%
%   INPUT:
%
%   N=total number of sample points
%   rL=rough surface length
%   h=rms height
%   kl=lower wavenumber cutoff
%   ku=upper wavenumber cutoff
%   Nf=number of tones
%   sdim=fractal dimension (1 <= s < 2)
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
Phi=2*pi*randn(Nf,1);
b=(ku/kl)^(1/(Nf-1));
Cn=sqrt(2*(1-b^(2*(sdim-2)))/(1-b^(2*(sdim-2)*Nf)));

dx=rL/N;
x=[-N/2+1:1:N/2]*dx;

f=sin(kl*x+Phi(1));
for n=1:Nf-1
	f=f+b^((sdim-2)*n)*sin(kl*b^n*x+Phi(n+1));
end
f=f*h*Cn;

n=2:N-1;
df1=(f(n+1)-f(n-1))/(2*dx);
df=[(f(2)-f(N))/(2*dx),df1,(f(1)-f(N-1))/(2*dx)];
