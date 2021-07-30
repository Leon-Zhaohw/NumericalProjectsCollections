%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             nsphere.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scattering by N Small Spheres, dipole point interactions only
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: By Chite Chen, April 1998

clear all;

% Input Parameters
Nr=30;            % Define Number of realization (Nr)
N=200;            % Define Number of particles (N)
f=0.1;            % Define fractional volume
L=1;              % Define the size of box
ka=0.2;           % Given ka value
eps_p=3.2;        % Define relative permittivity of particle 

x_p=eps_p-1;

% Define the step observation angle theta
step=2;           % in degree

% Calculate the particle size
a=((3*f*L^3)/(4*pi*N))^(1/3);
v=(4*pi/3*a^3);

% Calculate wavenumber and wavelength
k=ka/a;
lambda=2*pi/k;

load pos.dat;       % from mcpdf.m

for i_r=1:Nr,       % index for realization;
  location=position((i_r-1)*N+1:(i_r*N),:);
  x=location(:,1);
  y=location(:,2);
  z=location(:,3);
  X=x*ones(1,N);
  Y=y*ones(1,N);
  Z=z*ones(1,N);
  R=sqrt((X-X.').^2+(Y-Y.').^2+(Z-Z.').^2);
%
% let the diagonal of R be unzero (K.H. Ding)
  R=R+eye(N);       
%
  G1=(-1+i*k*R+k^2*R.^2).*exp(i*k*R)./(4*pi*k^2*R.^3);
  G2=(3-3*i*k*R-k^2*R.^2).*exp(i*k*R)./(4*pi*k^2*R.^5);

% Solving for [a] which [B][A] = [V]
% Matrix filling for [B]
  Bxx=-k^2*x_p*v*[G1+(X-X.').^2.*G2];
  Bxy=-k^2*x_p*v*((X-X.').*(Y-Y.').*G2);
  Bxz=-k^2*x_p*v*((X-X.').*(Z-Z.').*G2);
  Byx=Bxy;
  Byy=-k^2*x_p*v*[G1+(Y-Y.').^2.*G2];
  Byz=-k^2*x_p*v*((Y-Y.').*(Z-Z.').*G2);
  Bzx=Bxz;
  Bzy=Byz;
  Bzz=-k^2*x_p*v*[G1+(Z-Z.').^2.*G2];

% For self patch terms (diagonal elements in each B matrix)
  for mm=1:N,
    Bxx(mm,mm)=1+x_p/3;
    Bxy(mm,mm)=0;
    Bxz(mm,mm)=0;
    Byx(mm,mm)=0;
    Byy(mm,mm)=1+x_p/3;
    Byz(mm,mm)=0;
    Bzx(mm,mm)=0;
    Bzy(mm,mm)=0;
    Bzz(mm,mm)=1+x_p/3;
  end
  B=[Bxx Bxy Bxz; Byx Byy Byz; Bzx Bzy Bzz];

% Fill in [V]
  V=zeros(3*N,1);
  V(1:N)=sqrt(v)*exp(i*k*z);
  A=B\V;

% Calculating Scattering Field
  theta=0:step:180;          % observation angle in degree
  phi=0;
  theta=theta*pi/180;        % observation angle in radian
  phi=phi*pi/180;
  n_theta=size(theta,2);
  for i_theta=1: n_theta;
    ks=[sin(theta(i_theta))*cos(phi) sin(theta(i_theta))*sin(phi) cos(theta(i_theta))];
    vs=[cos(theta(i_theta))*cos(phi) cos(theta(i_theta))*sin(phi) -sin(theta(i_theta))];
    hs=[-sin(phi) cos(phi) 0];
    phase=exp(-i*location*k*ks.');
    C_Evs=[vs(1)*phase.' vs(2)*phase.' vs(3)*phase.'];
    Evs(i_r,i_theta)=k^2/(4*pi)*x_p*sqrt(v)*(C_Evs*A);
    C_Ehs=[hs(1)*phase.' hs(2)*phase.' hs(3)*phase.'];
    Ehs(i_r,i_theta)=k^2/(4*pi)*x_p*sqrt(v)*(C_Ehs*A);
  end
end

% Calculate <abs(incoherent filed)^2>
phase_fun_v=sum((abs(Evs)).^2)/Nr-(abs(sum(Evs)/Nr)).^2;
phase_fun_h=sum((abs(Ehs)).^2)/Nr-(abs(sum(Ehs)/Nr)).^2;

figure;
subplot(211),plot(theta*180/pi,phase_fun_v);
grid on;
title('Phase Function Evs; ka=0.2, phi=0')
legend('Nr=30, N=200')
xlabel('\theta in degree');
subplot(212),plot(theta*180/pi,phase_fun_h);
grid on;
title('Phase Function Ehs')
legend('Nr=30, N=200')
xlabel('\theta in degree');