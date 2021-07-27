%
%   Convect a Gaussian pulse using plane rotation on [-1,1]^2.
%
%   Usage:   
%
%   N=33; c_acc
%
%   This version has an O(dt^2) accurate initial step, but with
%   a very small initial dt on the first revolution (which consequently
%   requires a couple of extra steps.
%
%   See cfast.m for a simpler startup phase, which is normally suitable
%   for most fluid simulations.
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

CFL = .500;                                      % CFL constraint ( < ~.5 )
dxmin = min(abs(diff(x))); Dt = CFL*dxmin;    
Nsteps = ceil(Tfinal/Dt);  Dt = Tfinal/Nsteps;   % Integer # steps/period
[Dt Nsteps]

u=u0; u1=0*u; u2=u1; e1=u1; e2=u1;   % Initialize lagged fields

Npass=20; time=0; 

for pass=1:Npass;

  istart=1; istop=Nsteps;
  if pass==1;  istart=-1; end; % SPECIAL START UP

  for step=istart:Nsteps; 

    if pass==1 && step==istart;
       dt = .01*Dt; t_bdf = [time];
    elseif pass==1 && step==istart+1;
       dt = .10*Dt;
    elseif pass==1 && step==istart+2;
       dt = Dt-time;
    else
       dt = Dt;
    end;

    time=time+dt;
    t_bdf = [time t_bdf]; kmax=3; k=length(t_bdf)-1; k=min(k,kmax);
    t_bdf = t_bdf(1:k+1); t_ext = t_bdf(2:k+1);
    a=zeros(kmax+1,1); b=zeros(kmax+1,1);

    aw=fd_weights_full(time,t_ext,0); a(1:k) = aw(1:k); % Interpolant at t=time
    bw=fd_weights_full(time,t_bdf,1); b(1:k)=bw(2:k+1,2);b0=bw(1,2); % Deriv @ t

    ux = cxo.*(Dx*u); uy = cyo.*(u*Dy');
    e3=e2; e2=e1; e1=ux+uy;

    u3=u2; u2=u1; u1=u;

    u = -( u1*b(1)+u2*b(2)+u3*b(3) + e3*a(3)+e2*a(2)+e1*a(1) )/b0;

%   if mod(step,1000)==0;mesh(xx,yy,u); axis([-1 1 -1 1 -.1 1]); pause(.1); end;
%   if mod(step,1)==0;mesh(xx,yy,u); axis([-1 1 -1 1 -.1 1]); pause; end;

  end;

  mesh(xx,yy,u); axis([-1 1 -1 1 -.1 1]); pause(.1);
  err(pass) = max(max(abs(u-u0))); t(pass) = time;

end;

figure
semilogy(t,err,'ro',t,err,'r-'); axis square;
xlabel('time'); ylabel('maximum pointwise error');
title('Error for Convection, Acc. Start')


