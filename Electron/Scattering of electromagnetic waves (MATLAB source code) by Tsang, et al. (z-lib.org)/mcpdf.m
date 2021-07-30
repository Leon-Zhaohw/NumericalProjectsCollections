%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             mcpdf.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main program for the Monte Carlo simulations of pair distribution 
% functions and the creation of a series of random realizations of 
% particle positions
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: by K.H. Ding, November 1998

% Input Parameters
ntot=input('Enter total number of spheres : ');
fv=input('Enter fractional volume of spheres (< 0.4) : ');
cnst=input('Enter maximum displacement (< 1) : ');
npsr=input('Enter number of passes for each realization : ');
nrlz=input('Enter total number of realizations : ');
seed=input('Enter seed for random numbers : ');

vol=1.0;
rho=ntot/vol;
da=(6.0*fv/pi/rho)^(1.0/3.0);
nd=fix(1.0/da);
ncell=nd*nd*nd;
dl=1.0/nd;
ntpas=npsr*nrlz;
dinc=1.0/ntot/ntpas;
rgmax=5.;
srgmax=rgmax*da;
if srgmax >= 0.5
  srgmax=0.5;
end
rgmax=srgmax/da;
del=cnst*da;

if ncell < ntot
  fprintf('\n Number of Spheres > Number of Cells ==> STOP ! \n');
  break
end

if dl <= da
  fprintf('\n Diameter > Cell Length ==> STOP ! \n');
  break
end

fpos=fopen('pos.dat','w+');
fpdf=fopen('pdf.dat','w+');

% initial regular setting of spheres
xrow=zeros(1,ntot);
yrow=zeros(1,ntot);
zrow=zeros(1,ntot);
da2=ones(ntot,ntot);
da2=da*da*da2;
np=0;
for i=0:nd-1
  if np > ntot
    break
  else
    for j=0:nd-1
      if np > ntot
        break
      else
        for k=0:nd-1
          np=np+1;
          if np > ntot
            break
          else
            xrow(np)=(k+0.5)*dl;
            yrow(np)=(j+0.5)*dl;
            zrow(np)=(i+0.5)*dl;
          end
        end
      end
    end
  end
end

mg=40;
ddcst=4.0*pi*(ntot-1.0)*rho*da*da*da/3.0;
dgr=(rgmax-1.0)/mg;
r=zeros(1,mg);
ddr=zeros(1,mg);
f=zeros(1,mg);
fmn=zeros(1,mg);
g=zeros(1,mg);
for n=1:mg
  rl=1.0+n*dgr;
  rs=rl-dgr;
  r(n)=sqrt((rl*rl+rs*rs)/2.0);
  ddr(n)=rl*rl*rl-rs*rs*rs;
end

rand('seed',seed);

erow=ones(1,ntot);
ecol=ones(ntot,1);
move=zeros(ntot,1);
stay=ones(ntot,1);
rx=zeros(ntot,ntot);
ry=zeros(ntot,ntot);
rz=zeros(ntot,ntot);
rr2=zeros(ntot,ntot);
ix=zeros(ntot,ntot);

% initial random setting shuffling
fprintf('\n Initial Shuffling Starts ... \n');
acp=0.0;
ovp=0.0;
for ir=1:npsr
  ranx=erow-2.0*rand(1,ntot);
  rany=erow-2.0*rand(1,ntot);
  ranz=erow-2.0*rand(1,ntot);
  xtry=xrow+del*ranx;
  ytry=yrow+del*rany;
  ztry=zrow+del*ranz;
  xtry=xtry-fix(2.0*xtry-erow);
  ytry=ytry-fix(2.0*ytry-erow);
  ztry=ztry-fix(2.0*ztry-erow);
  xtry=xtry-fix(2.0*xtry-erow);
  ytry=ytry-fix(2.0*ytry-erow);
  ztry=ztry-fix(2.0*ztry-erow);

% check separation between pairs of spheres
  rx=xtry'*erow-ecol*xrow;
  ry=ytry'*erow-ecol*yrow;
  rz=ztry'*erow-ecol*zrow;
  rx=rx-fix(2.0*rx);
  ry=ry-fix(2.0*ry);
  rz=rz-fix(2.0*rz);
  rr2=rx.^2+ry.^2+rz.^2;
  for i=1:ntot
     rr2(i,i)=da*da;
  end
  gt=rr2 >= da2;
  for i=1:ntot
    move(i)=all(gt(i,:));
  end
  stay=ecol-move;
  xrow=move'.*xtry+stay'.*xrow;
  yrow=move'.*ytry+stay'.*yrow;
  zrow=move'.*ztry+stay'.*zrow;
  for i=1:ntot
    acp=acp+dinc*move(i);
    ovp=ovp+dinc*stay(i);
  end
end
fprintf('\n acceptance rate = %8.4f \n',acp);
fprintf(' overlaping rate = %8.4f \n',ovp);
fprintf('\n Initial Shuffling Done !!! \n');

% Monte Carlo shuffling
ir=0;
acp=0.0;
ovp=0.0;

fprintf('\n Monte Carlo Shuffling Starts ... \n');
for ip=1:ntpas
  ipp=ip-npsr;
  while ipp > 0
    ipp=ipp-npsr;
  end

  ranx=erow-2.0*rand(1,ntot);
  rany=erow-2.0*rand(1,ntot);
  ranz=erow-2.0*rand(1,ntot);
  xtry=xrow+del*ranx;
  ytry=yrow+del*rany;
  ztry=zrow+del*ranz;
  xtry=xtry-fix(2.0*xtry-erow);
  ytry=ytry-fix(2.0*ytry-erow);
  ztry=ztry-fix(2.0*ztry-erow);
  xtry=xtry-fix(2.0*xtry-erow);
  ytry=ytry-fix(2.0*ytry-erow);
  ztry=ztry-fix(2.0*ztry-erow);

% check separation between pairs of spheres
  rx=xtry'*erow-ecol*xrow;
  ry=ytry'*erow-ecol*yrow;
  rz=ztry'*erow-ecol*zrow;
  rx=rx-fix(2.0*rx);
  ry=ry-fix(2.0*ry);
  rz=rz-fix(2.0*rz);
  rr2=rx.^2+ry.^2+rz.^2;
  for i=1:ntot
    rr2(i,i)=da*da;
  end
  gt=rr2 >= da2;
  for i=1:ntot
    move(i)=all(gt(i,:));
  end
  stay=ecol-move;
  xrow=move'.*xtry+stay'.*xrow;
  yrow=move'.*ytry+stay'.*yrow;
  zrow=move'.*ztry+stay'.*zrow;
  for i=1:ntot
    acp=acp+dinc*move(i);
    ovp=ovp+dinc*stay(i);
  end

% tabulate the occurance of pair separations
  fmn=f;
  rx=xrow'*erow-ecol*xrow;
  ry=yrow'*erow-ecol*yrow;
  rz=zrow'*erow-ecol*zrow;
  rx=rx-fix(2.0*rx);
  ry=ry-fix(2.0*ry);
  rz=rz-fix(2.0*rz);
  rr2=rx.^2+ry.^2+rz.^2;
  ix=fix((sqrt(rr2/da/da)-1.0)/dgr)+1;
  for i=1:ntot
    for j=1:ntot
      if ix(i,j) >= 1 & ix(i,j) <= mg
        fmn(ix(i,j))=fmn(ix(i,j))+1.0;
      end
    end
  end
  f=fmn;

  if ipp == 0
    ir=ir+1;
    fprintf('\n     realization = %6u \n',ir);
    fprintf('            pass = %6u \n',ip);
    fprintf(' acceptance rate = %8.4f \n',acp);
    fprintf(' overlaping rate = %8.4f \n',ovp);
%   output position
    fprintf(fpos,'%6u \n',ir);
    for np=1:ntot
      fprintf(fpos,'%8.4f %8.4f %8.4f \n',xrow(np),yrow(np),zrow(np));
    end
  end
end

fprintf(' \n\n');
fprintf(' total acceptance rate = %8.4f \n',acp);
fprintf(' total overlaping rate = %8.4f \n',ovp);
fprintf('\n Monte Carlo Shuffling Done !!! \n');

% output pair distribution function
for jj=1:mg
  g(jj)=f(jj)/ddcst/ddr(jj)/ntpas;
  fprintf(fpdf,'%8.4f %8.4f \n',r(jj),g(jj));
end
fclose(fpos);
fclose(fpdf);

% plot pair distribution function
% (requires output of pypdf.m)
% clear;
% load pdf.dat -ascii;
% figure;
% plot(pdf(:,1),pdf(:,2),'ro');
% axis([0.0,5.0,0.0,4.0]);
% xlabel('r/d ');
% ylabel('g(r) ');
% legend('Monte Carlo');