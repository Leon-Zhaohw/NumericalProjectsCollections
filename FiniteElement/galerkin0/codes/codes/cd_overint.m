
%
%   Convect a Gaussian pulse using plane rotation on [-1,1]^2 with:
%
%   .stable formulation with overintegration for convection term (dealiasing)
%   .no filtering
%   .diffusion
%   .all Dirichlet conditions
%   .single domain
%
%
%   Usage:   
%
%   N=33; cd_overint
%
%   The purpose here is to demonstrate a stable convection/diffusion formulation
%
%   This version has only an O(dt^2) accurate initial step, but
%   this appears to be adequate in most cases.
%

%%%
%%%   INITIALIZATION
%%%

[Ah,Bh,Ch,Dh,z,w] = SEMhat(N); 

No = floor(3*N/2);                 % Overintgration for Convection
[Ao,Bo,Co,Do,zo,w] = SEMhat(No); 
J  = interp_mat(zo,z);             % Interpolator to fine mesh
bo = diag(Bo); Bo=bo*bo';          % Diagonal mass matrix on fine mesh

x=z; Dx=Dh; nx=size(Ah,1); 
y=z; Dy=Dh; ny=size(Ah,1); 

%  Set up interpolation matrix to uniform mesh for plotting purposes
nf=ceil(1.5*N); zf=(0:nf)./nf; zf= -1 + 2*zf; 
Jf=Interp_mat(zf,z);
[xf,yf]=meshgrid(zf,zf); xf=xf'; yf=yf';

nx=size(Ah,1); x=z; 
Ax=eye(nx); Ax(2:nx-1,2:nx-1)=Ah(2:nx-1,2:nx-1);  % Interior only, for
Bx=eye(nx); Bx(2:nx-1,2:nx-1)=Bh(2:nx-1,2:nx-1);  % Dirichlet conditions
bx=diag(Bx); Bx=sparse(Bx); Dx=Dh;

y=x; ny=nx; Ay=Ax; By=Bx; Dy=Dx; by=bx;


%  Define Mesh
[xx,yy] = meshgrid(x,y); xx=xx';yy=yy'; % Lexicographical ordering

%  Diagonal mask (1 in interior, 0 on boundary) to enforce u=0 on dOmega.
mx=ones(nx,1); mx(1)=0; mx(nx)=0;
my=ones(ny,1); my(1)=0; my(ny)=0;
mask = mx*my';

%  Diagonal mass matrix, with mask included.
mass = mask.*(bx*by');

%
%  Diagonalize 1D Poisson operators
%
Ax=full(Ax); Bx=full(Bx); [Sx,lx] = eig(Ax,Bx); 
Ay=full(Ay); By=full(By); [Sy,ly] = eig(Ay,By); 
   
%  Normalize columns of S  ( MATLAB is _not_ LAPACK )
scalex = Sx'*Bx*Sx; scalex=diag(1./sqrt(diag(scalex))); Sx = Sx*scalex;
scaley = Sy'*By*Sy; scaley=diag(1./sqrt(diag(scaley))); Sy = Sy*scaley;
   
   
%  Eigenvalues for fast solver:  I x Lambda + Lambda x I
ix = ones(nx,1);    lx = diag(lx);      
iy = ones(ny,1);    ly = diag(ly);


%  Define velocity  (comment/uncomment to test other velocities)
pi2=pi/2;
 cx = -cos(pi2*xx).*sin(pi2*yy); cy =  sin(pi2*xx).*cos(pi2*yy); 
%cx = .5+0*yy;                   cy=1-0*xx; 
%cx = -yy;                       cy=xx; 

%  Velocity on fine mesh, collocated w/ overintegration weights
cxo = Bo.*(J*cx*J'); cyo = Bo.*(J*cy*J');


%  Set initial condition; enforce zero on boundary w/ quadratic bubble
x0 =-0.0; y0=-0.5; delta = 0.10; rr=(xx-x0).^2+(yy-y0).^2;
u0 = exp(-((rr./(delta^2)).^1)).*(1-xx.*xx).*(1-yy.*yy);
mesh(xx,yy,u0); axis([-1 1 -1 1 -.1 1]); pause(.1);

%  Define viscosity, Tfinal, & CFL

nu = .0001;                                      % Viscosity
Tfinal = 2*pi;                                   % Time for single revolution

CFL = .5;                                        % CFL constraint ( < ~.5 )
dxmin = min(abs(diff(x))); Dt = CFL*dxmin;    
Nsteps = ceil(Tfinal/Dt);  Dt = Tfinal/Nsteps;   % Integer # steps/period
format long; [Dt Nsteps]


Npass=20; time=0; dt=Dt;

u=u0; u1=0*u; u2=u1; e1=u1; e2=u1;   % Initialize lagged fields

%%%
%%%   MAIN TIMESTEPPING LOOP
%%%

for pass=1:Npass;  
  for step=1:Nsteps;

    if     pass==1 && step==1; % Fixed dt, only O(dt^2) accuracy at start
       b0=1/dt; b=([ -1 0 0 ]')/dt;                 a=[ 1  0 0 ]'; 
       dxy = nu*(lx*iy' + ix*ly') + b0*(ix*iy'); dxyi = 1./dxy; % H = nu*A + b0*B
    elseif pass==1 && step==2; 
       b0=1.5/dt; b=([ -4 1 0 ]')./(2*dt);          a=[ 2 -1 0 ]';
       dxy = nu*(lx*iy' + ix*ly') + b0*(ix*iy'); dxyi = 1./dxy;
    elseif pass==1 && step==3; 
       b0=11./(6*dt); b=[ -18 9 -2 ]; b=b'./(6*dt); a=[ 3 -3 1 ]';
       dxy = nu*(lx*iy' + ix*ly') + b0*(ix*iy'); dxyi = 1./dxy;
    end;

    time = time+dt;

    ux = cxo.*(J*Dx*u*J'); uy = cyo.*(J*u*Dy'*J'); % Mass matrix included
    e3=e2; e2=e1; e1 = J'*(ux+uy)*J;

    u3=u2; u2=u1; u1=u;

%   Mask included in mass matrix
    r = -mass.*(u1*b(1)+u2*b(2)+u3*b(3))-mask.*(e1*a(1)+e2*a(2)+e3*a(3));
    u = Sx*( (Sx'*r*Sy).*dxyi )*Sy';   %  u = Hinv f


%   if mod(step,200)==0;mesh(xx,yy,u);axis([-1 1 -1 1 -.1 1]);pause(.1);end;
%   if mod(step,200)==0;contour(xx,yy,u,25);axis square; pause(.1);end;
    if mod(step,200)==0; uf=Jf*u*Jf';  % Interpolate onto fine uniform mesh for contour plot
      contour(xf,yf,uf,25);axis square; pause(.1);
    end;

  end;

end;
