%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              pypdf.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main program for the calculation of Percus Yevick pair distribution 
% function and structure factor for a medium with spherical particles 
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: by K.H. Ding, November 1998

% input parameters
fv=input('Enter the volume fraction of particles : ');
dia=input('Enter the diameter of particles : ');

tau= 999999999.
tauc=(2-sqrt(2))/6;
xic=(3*sqrt(2)-4)/2;
rm=5;
nk=512;
dk=pi/(rm*dia);
dr=rm*dia/(nk+1);
hk=zeros(2*nk+2,1);

vol=pi*dia^3/6;
rho=fv/vol;
xi=fv;
xi1=1-xi;
nu=tau+xi/xi1;
gama=xi*(1+xi/2)/(3*xi1^2);
lamd=6*(nu-sqrt(nu^2-gama))/xi;
mu=lamd*xi*xi1;
if mu > 1+2*xi
  break;
end

cst1=xi/xi1;
cst2=1-lamd*xi+3*cst1;
cst3=3-lamd*xi1;

fpk=fopen('pysf.dat','w+');

pk=0;
pp1=cst1*(4-lamd+3*cst1)+1;
pp2=0;
pyhk=(1/(pp1^2+pp2^2)-1)/rho;
pysf=1+pyhk*rho;
fprintf(fpk,'%6u %14.9f \n',pk,pysf);

for ik=1:nk
  pk=ik*dk;
  x=pk*dia/2;
  snx=sin(x);
  csx=cos(x);
  psix=snx/x;
  phix=3*(snx-x*csx)/(x^3);
  pp1=cst1*(cst2*phix+cst3*psix)+csx;
  pp2=cst1*x*phix+snx;
  pyhk=(1/(pp1^2+pp2^2)-1)/rho;
  pysf=1+pyhk*rho;
  fprintf(fpk,'%6u %14.9f \n',pk,pysf);
  hk(ik+1)=pk*pyhk;
end
fclose(fpk);

hw=-2*fft(hk);
hr=imag(hw(2:nk+1));

fpr=fopen('pypdf.dat','w+');
for ir=1:nk
  r=ir*dr/dia;
  if r >= 1
    g=1+hr(ir)/(rm*4*pi*r*dia^2);
    fprintf(fpr,'%14.9f %14.9f \n',r,g);
  end
end
fclose(fpr);

% plot pair distribution function
clear;
load pdfpy.dat -ascii;
load pysf.dat -ascii;
figure;
plot(pypdf(:,1),pypdf(:,2),'r:');
axis([0.0,5.0,0.0,3.0]);
title('PY pair distribution function');
xlabel('r/(2a) ');
ylabel('g(r) ');
legend('fv=0.3');

figure;
plot(pysf(:,1),pysf(:,2),'r:');
axis([0.0,150.0,0.0,2.0]);
title('PY structure factor');
xlabel('k ');
ylabel('S(k) ');
legend('fv=0.3');
clear;
