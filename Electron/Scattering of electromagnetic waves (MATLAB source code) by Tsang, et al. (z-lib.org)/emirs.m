function epol=emirs(ipol,theta,lambda,epsilon,lc,hh)
%EMIRS computes the emissivity from 2-D Gaussian dielectric rough 
%surface using small perturbation method (SPM)
%
%   epol=EMIRS(ipol,theta,lambda,epsilon,lc,hh)
%
%   INPUT:
%
%   ipol=calculate TE (1) or TM (2)emission
%   theta(npoints)=observed polar angles (deg) (0 to 90)
%   lambda=wavelength in free space
%   epsilon=relative permittivity of lower medium
%   lc=rough surface correlation length in lambda
%   hh=rough surface rms height in lambda
%
%   OUTPUT:
%
%   epol(npoints)=total emissivity in TE or TM for each angle theta
%
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: by Chite Chen, November 1998

k=2*pi/lambda;
k1=k*sqrt(epsilon);

l=lc*lambda;
h=hh*lambda;
npoints=length(theta);

% currently only consider phi=0
phi=zeros(size(theta));

theta_i=theta*pi/180;
phi_i=phi*pi/180; 
kxi=k*sin(theta_i)*cos(phi_i);
kyi=k*sin(theta_i)*sin(phi_i);
kzi=k*cos(theta_i);
k1zi=sqrt(k1^2-(k*sin(theta_i)).^2);

% TE
if (ipol == 1) 
  term1=abs(Rh0(theta_i,k,k1)).^2;
  for ii=1:npoints,
    thetai=theta_i(ii);
    phii=phi_i(ii);
    term2(ii)=2*real(Rh0(thetai,k,k1)*conj(f_ee_2(thetai,phii,k,k1,h,l)));
    term3(ii)=eh3(thetai,phii,k,k1,h,l);
    term4(ii)=eh4(thetai,phii,k,k1,h,l);
  end
end

% TM
if (ipol == 2) 
  term1=abs(Rv0(theta_i,k,k1)).^2;
  for ii=1:npoints,
    thetai=theta_i(ii);
    phii=phi_i(ii);
    term2(ii)=2*real(Rv0(thetai,k,k1)*conj(f_hh_2(thetai,phi_i,k,k1,h,l)));
    term3(ii)=ev3(thetai,phii,k,k1,h,l);
    term4(ii)=ev4(thetai,phii,k,k1,h,l);
  end
end

epol=1-term1-term2-term3-term4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           W_Spectrum.m                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ww]=W_spectrum(kx,ky,k,k1,h,l)
% define the spectrum of interest (Gaussian spectral density)

ww=(h*l)^2/(4*pi)*exp(-((kx.^2+ky.^2)/4*l^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Rh0.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rh0]=Rh0(theta_i,k,k1)
% 0th order solution: reflection coefficient for h-polarization

kzi=k*cos(theta_i);
k1zi=sqrt(k1^2-(k*sin(theta_i)).^2);
Rh0=(kzi-k1zi)./(kzi+k1zi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Rv0.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rv0]=Rv0(theta_i,k,k1)
% 0th order solution: reflection coefficient for v-polarization

kzi=k*cos(theta_i);
k1zi=sqrt(k1^2-(k*sin(theta_i)).^2);
Rv0=(k1^2*kzi-k^2*k1zi)./(k1^2*kzi+k^2*k1zi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_hh_1.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=f_hh_1(thetai,phi_i,thetak,phi_k,k,k1,h,l)
% calculate f_hh_1 for a given spectrum

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
krhoi=k*sin(thetai);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);
kz=k*cos(thetak);
krho=k*sin(thetak);
k1z=sqrt(k1^2-(k*sin(thetak))^2);
y1=-k1z*k1zi.*cos(phi_k-phi_i);
y2=krho*krhoi*k1^2/k^2;
output=(k1^2-k^2)./(k1^2*kz+k^2*k1z)*2*k^2*kzi/(k1^2*kzi+k^2*k1zi).*(y1+y2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_he_1.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=f_he_1(thetai,phi_i,thetak,phi_k,k,k1,h,l)
% calculate f_he_1 for a given spectrum

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
krhoi=k*sin(thetai);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);
kz=k*cos(thetak);
krho=k*sin(thetak);
k1z=sqrt(k1^2-(k*sin(thetak))^2);
output=(k1^2-k^2)*k1z*k./(k1^2*kz+k^2*k1z)*2*kzi/(kzi+k1zi).*sin(phi_k-phi_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_eh_1.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=f_eh_1(thetai,phi_i,thetak,phi_k,k,k1,h,l)
% calculate f_eh_1 for a given spectrum

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);
kz=k*cos(thetak);
k1z=sqrt(k1^2-(k*sin(thetak))^2);
output=(k1^2-k^2)./(kz+k1z)*kzi*2*k*k1zi/(k^2*k1zi+k1^2*kzi).*sin(phi_k-phi_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_ee_1.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=f_ee_1(thetai,phi_i,thetak,phi_k,k,k1,h,l)
% calculate f_ee_1 for a given spectrum

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);
kz=k*cos(thetak);
k1z=sqrt(k1^2-(k*sin(thetak))^2);
output=(k1^2-k^2)./(kz+k1z)*2*kzi/(kzi+k1zi).*cos(phi_k-phi_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_hh_2.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=f_hh_2(thetai,phi_i,k,k1,h,l)
% calculate f_hh_2 for a given spectrum

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
krhoi=k*sin(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

% integrate the first part of f_hh_2 using numerical integration
n_kx=400;
n_ky=400;
aa=5/l;

kx=linspace(kxi-aa,kxi+aa,n_kx);
ky=linspace(kyi-aa,kyi+aa,n_ky);
delta_kx=2*aa/n_kx;
delta_ky=2*aa/n_ky;

for iky=1:n_ky,
  kyy=ky(iky);
  ww=W_spectrum(kx-kxi,kyy-kyi,k,k1,h,l);
  y11(iky)=trapz(kx,ww);
end
y11_sum = trapz(ky,y11);
y1_sum = k1^2*k*kzi^2*y11_sum;

% integrate the second part of f_hh_2 using numerical integration
for iky=1:n_ky,
  kyy=ky(iky);
  phi_k=atan2(kyy,kx);              % arc tangent
  krho=sqrt(kx.^2+kyy^2);
  kz=sqrt(k^2-krho.^2);
  k1z=sqrt(k1^2-krho.^2);
  ww=W_spectrum(kx-kxi,kyy-kyi,k,k1,h,l);
  y21=-k*k1zi^2*(k1^2-k^2)./(k1z+kz).*(sin(phi_k-phi_i)).^2;
  y22=k1^2/k*(k1^2-k^2)*krhoi^2*krho.^2;
  y23=2*krhoi*krho*k*k1^2.*(k1z+kz)*k1zi.*cos(phi_k-phi_i);
  y24=k*kz.*k1z*(k1^2-k^2)*k1zi^2.*(cos(phi_k-phi_i)).^2;
  y25=1./(k1^2*kz+k^2*k1z).*(y22-y23-y24);
  y_ww=(y21+y25).*ww;
  y2(iky)=trapz(kx,y_ww);
end
y22_sum=trapz(ky,y2);
y2_sum=y22_sum;
y_sum=k1zi*y1_sum+kzi^2*y2_sum;
output=-2*k*(k1^2-k^2)/(kzi*(k1^2*kzi+k^2*k1zi).^2)*(y_sum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_ee_2.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=f_ee_2(thetai,phi_i,k,k1,h,l)
% calculate f_ee_2 for a given spectrum

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
krhoi=k*sin(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

% integrate the first part of f_ee_2 using numerical integration
n_kx=400;
n_ky=400;
aa=5/l;

kx=linspace(kxi-aa,kxi+aa,n_kx);
ky=linspace(kyi-aa,kyi+aa,n_ky);
delta_kx=2*aa/n_kx;
delta_ky=2*aa/n_ky;

for iky=1:n_ky,
  kyy=ky(iky);
  ww=W_spectrum(kx-kxi,kyy-kyi,k,k1,h,l);
  y11(iky)=trapz(kx,ww);
end
y11_sum=trapz(ky,y11);
y1_sum=y11_sum;

% integrate the second part of f_ee_2 using numerical integration
for iky=1:n_ky,
  kyy=ky(iky);
  phi_k=atan2(kyy,kx);              % arc tangent
  krho=sqrt(kx.^2+kyy^2);
  kz=sqrt(k^2-krho.^2);
  k1z=sqrt(k1^2-krho.^2);
  ww=W_spectrum(kx-kxi,kyy-kyi,k,k1,h,l);
  y21=kz.*k1z./(k1^2*kz+k^2*k1z).*(sin(phi_k-phi_i)).^2;
  y22=(cos(phi_k-phi_i)).^2./(k1z+kz);
  y_ww=(y21+y22).*ww;
  y2(iky)=trapz(kx,y_ww);
end
y22_sum=trapz(ky,y2);
y2_sum=y22_sum;
y_sum=-k1zi*y1_sum+(k1^2-k^2)*y2_sum;
output=-2*kzi*(k1^2-k^2)/(kzi+k1zi)^2*(y_sum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               eh3.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=eh3(thetai,phi_i,k,k1,h,l)
% calculate the integration of f_ee_1 using trapezoidal rule

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

n_theta_k_points=180;
n_phi_k_points=720;

theta_k=linspace(0,pi/2,n_theta_k_points);
phi_k=linspace(0,2*pi,n_phi_k_points);
delta_theta_k=(pi/2-0)/n_theta_k_points;
delta_phi_k=(2*pi-0)/n_phi_k_points;

for ik=1:n_theta_k_points,
  thetak=theta_k(ik);
  kx=k*sin(thetak)*cos(phi_k);
  ky=k*sin(thetak)*sin(phi_k);
  kz=k*cos(thetak);
  k1z=sqrt(k1^2-(k*sin(thetak))^2);
  ww=W_spectrum(kx-kxi,ky-kyi,k,k1,h,l);
  fee=f_ee_1(thetai,phi_i,thetak,phi_k,k,k1,h,l);
  y=ww.*(abs(fee)).^2;
  y_theta_k_sum=trapz(phi_k,y);
% the inner integral value for a given theta_k
  y_sum(ik)=sin(thetak)*(cos(thetak)).^2*y_theta_k_sum;
end
yy=trapz(theta_k,y_sum);
output=k^2/cos(thetai)*yy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               eh4.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=eh4(thetai,phi_i,k,k1,h,l)
% calculate the integration of f_he_1 using trapezoidal rule

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

n_theta_k_points=180;
n_phi_k_points=720;

theta_k=linspace(0,pi/2,n_theta_k_points);
phi_k=linspace(0,2*pi,n_phi_k_points);
delta_theta_k=(pi/2-0)/n_theta_k_points;
delta_phi_k=(2*pi-0)/n_phi_k_points;

for ik=1:n_theta_k_points,
  thetak=theta_k(ik);
  kx=k*sin(thetak)*cos(phi_k);
  ky=k*sin(thetak)*sin(phi_k);
  kz=k*cos(thetak);
  k1z=sqrt(k1^2-(k*sin(thetak))^2);
  ww=W_spectrum(kx-kxi,ky-kyi,k,k1,h,l);
  fhe=f_he_1(thetai,phi_i,thetak,phi_k,k,k1,h,l);
  y=ww.*(abs(fhe)).^2;
  y_theta_k_sum=trapz(phi_k,y);
% the inner integral value for a given theta_k
  y_sum(ik)=sin(thetak)*(cos(thetak)).^2*y_theta_k_sum;
end
yy=trapz(theta_k,y_sum);
output=k^2/cos(thetai)*yy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ev3.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=ev3(thetai,phi_i,k,k1,h,l)
% calculate the integration of f_eh_1 using trapezoidal rule

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

n_theta_k_points=180;
n_phi_k_points=720;
theta_k=linspace(0,pi/2,n_theta_k_points);
phi_k=linspace(0,2*pi,n_phi_k_points);
delta_theta_k=(pi/2-0)/n_theta_k_points;
delta_phi_k=(2*pi-0)/n_phi_k_points;

for ik=1:n_theta_k_points,
  thetak=theta_k(ik);
  kx=k*sin(thetak)*cos(phi_k);
  ky=k*sin(thetak)*sin(phi_k);
  kz=k*cos(thetak);
  k1z=sqrt(k1^2-(k*sin(thetak))^2);
  ww=W_spectrum(kx-kxi,ky-kyi,k,k1,h,l);
  feh=f_eh_1(thetai,phi_i,thetak,phi_k,k,k1,h,l);
  y=ww.*(abs(feh)).^2;
  y_theta_k_sum=trapz(phi_k,y);
% the inner integral value for a given theta_k
  y_sum(ik)=sin(thetak)*(cos(thetak)).^2*y_theta_k_sum;
end
yy=trapz(theta_k,y_sum);
output=k^2/cos(thetai)*yy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ev4.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=ev4(thetai,phi_i,k,k1,h,l)
% calculate the integration of f_hh_1 using trapezoidal rule

kxi=k*sin(thetai)*cos(phi_i);
kyi=k*sin(thetai)*sin(phi_i);
kzi=k*cos(thetai);
k1zi=sqrt(k1^2-(k*sin(thetai))^2);

n_theta_k_points=180;
n_phi_k_points=720;
theta_k=linspace(0,pi/2,n_theta_k_points);
phi_k=linspace(0,2*pi,n_phi_k_points);
delta_theta_k=(pi/2-0)/n_theta_k_points;
delta_phi_k=(2*pi-0)/n_phi_k_points;

for ik=1:n_theta_k_points,
  thetak=theta_k(ik);
  kx=k*sin(thetak)*cos(phi_k);
  ky=k*sin(thetak)*sin(phi_k);
  kz=k*cos(thetak);
  k1z=sqrt(k1^2-(k*sin(thetak))^2);
  ww=W_spectrum(kx-kxi,ky-kyi,k,k1,h,l);
  fhh=f_hh_1(thetai,phi_i,thetak,phi_k,k,k1,h,l);
  y=ww.*(abs(fhh)).^2;
  y_theta_k_sum=trapz(phi_k,y);
% the inner integral value for a given theta_k
  y_sum(ik)=sin(thetak)*(cos(thetak)).^2*y_theta_k_sum;
end
yy=trapz(theta_k,y_sum);
output=k^2/cos(thetai)*yy;
