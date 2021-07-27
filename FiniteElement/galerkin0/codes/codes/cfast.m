%
%   Convect a Gaussian pulse using plane rotation on [-1,1]^2.
%
%   Usage:   
%
%   N=33; cfast
%
%   Single domain, explicit convection, no dealiasing, no filter
%
%   NOTE: w/o dealiasing, this case is restricted to
%         constant velocities or to rotational flows.
%
%   The purpose here is to demonstrate a minimalist SEM
%   implementation of convection that is as fast as possible.
%
%   For this particular case (single element), no mass matrix is required.
%
%   This version has only an O(dt^2) accurate initial step.
%   See c_acc.m for a more accurate startup phase.
%
%

[Ah,Bh,Ch,Dh,z,w] = SEMhat(N); 

x=z; Dx=Dh; nx=size(Ah,1); 
y=z; Dy=Dh; ny=size(Ah,1); 


% Diagonal mask (1 in interior, 0 on boundary)
mx=ones(nx,1); mx(1)=0; mx(nx)=0;
my=ones(ny,1); my(1)=0; my(ny)=0;
mask = mx*my';


%  Define Mesh
[xx,yy] = meshgrid(x,y); xx=xx';yy=yy'; % Lexicographical ordering

pi2=pi/2;    %  Define Velocity  (comment/uncomment to test other velocities)
%cx = -cos(pi2*xx).*sin(pi2*yy); cy =  sin(pi2*xx).*cos(pi2*yy); 
%cx = .5+0*yy;                   cy=1-0*xx; 
 cx = -yy;                       cy=xx; 


% Premultiply velocity by mask to enforce Dirichlet condition
 cxo = mask.*cx; cyo = mask.*cy;


% Set initial condition; enforce zero on boundary w/ quadratic bubble
x0 =-0.0; y0=-0.5; delta = 0.10; rr=(xx-x0).^2+(yy-y0).^2;
u0 = exp(-((rr./(delta^2)).^1)).*(1-xx.*xx).*(1-yy.*yy);
mesh(xx,yy,u0); axis([-1 1 -1 1 -.1 1]); pause(.1);

Tfinal = 2*pi;                                   % Time for single revolution

CFL = .5;                                        % CFL constraint ( < ~.5 )
dxmin = min(abs(diff(x))); Dt = CFL*dxmin;    
Nsteps = ceil(Tfinal/Dt);  Dt = Tfinal/Nsteps;   % Integer # steps/period
format long; [Dt Nsteps]


Npass=20; time=0; dt=Dt;

u=u0; u1=0*u; u2=u1; e1=u1; e2=u1;   % Initialize lagged fields

for pass=1:Npass;                    % Simple scheme, only
  for step=1:Nsteps;                 %    O(dt^2) accuracy at start


    if pass==1 && step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end;
    if pass==1 && step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end;
    if pass==1 && step==3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end;

    time = time+dt;

    b0i = 1./b0;

    ux = cxo.*(Dx*u); uy = cyo.*(u*Dy');
    e3=e2; e2=e1; e1=ux+uy;

    u3=u2; u2=u1; u1=u;

    u = -b0i*( u1*b(1)+u2*b(2)+u3*b(3) + dt*(e3*a(3)+e2*a(2)+e1*a(1)));

  end;

  mesh(xx,yy,u); axis([-1 1 -1 1 -.1 1]); pause(.1);
  err(pass) = max(max(abs(u-u0))); t(pass) = time; 
  [pass time err(pass)]

end;

figure
semilogy(t,err,'ro',t,err,'r-'); axis square;
xlabel('time'); ylabel('maximum pointwise error');
title('Error for Convection, Std. Start')

