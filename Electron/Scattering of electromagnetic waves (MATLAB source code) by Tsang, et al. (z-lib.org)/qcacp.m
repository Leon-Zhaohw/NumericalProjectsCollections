%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             qcacp.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main program for the calculation of effective propagation constant 
% under the quasicrystalline approximation with coherent potentail 
% approximation (QCA-CP) for a medium with particle size distribution 
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: by K.H. Ding, November 1998

% input parameters
freq=input('Enter the frequency (GHz) : ');
ep0=input('Enter the relative permittivity of background : ');
eps=input('Enter the relative permittivity of scatterer : ');
alf=input('Enter the parameter "P" of size distribution : ');
gma=input('Enter the parameter "Q" of size distribution : ');
am=input('Enter the mode radius of size distribution (cm) : ');
na=input('Enter the number of discretized radius sizes : ');

wavl=30/freq;
k0=2*pi/wavl;
k=k0*sqrt(ep0);
ks=k0*sqrt(eps);
mu=((ks/k)^2-1)/3;

% determine constant "b" in size distribution
b=alf/gma/(am^gma);

pi6=pi/6;
alf3=alf+3;
alf6=alf+6;
af4g=(alf+4)/gma;
af6g=alf6/gma;
I=eye(na);
niter=10;

% determine the range of radius
ac=(af6g/b)^(1.0/gma);
eaf6g=exp(af6g);
aend=1.0;
np=1000;
da=ac/np;
amin=ac;
for ip=1:np
  a=ac-da*ip;
  ratio=((a/ac)^alf6)*eaf6g*exp(-b*a^gma);
  if ratio < 1.0e-3
    break;
  end
end
amin=a;
da=(aend-ac)/np;
amax=ac;
for ip=1:np
  a=ac+da*ip;
  ratio=((a/ac)^alf6)*eaf6g*exp(-b*a^gma);
  if ratio < 1.0e-3
    break;
  end
end
amax=a;

% determine discrete diameters "dia"
xa=zeros(1,na);
xb=zeros(1,na);
dia=zeros(1,na);
da=(amax-amin)/na;
for ia=1:na
  xa(ia)=amin+da*(ia-1);
  xb(ia)=xa(ia)+da;
  dia(ia)=xa(ia)+xb(ia);
end

fp=fopen('qcacp.dat','w+');
for nf=1:201
  ftot=0.002*(nf-1);
  if ftot == 0.0
    ner=0.0;
    edc=0.0;
  else

%   determine constant "c" in size distribution
    c=3*ftot*gma*b^af4g/4/pi/gamma(af4g);

%   determine discrete volume fractions "fv", and number
%   densities "rho", and contants xi0, xi1,xi2, and xi3
    fv=zeros(1,na);
    rho=zeros(1,na);
    xi0=0.0;
    xi1=0.0;
    xi2=0.0;
    xi3=0.0;
    c4p3=c*4*pi/3;
    for ia=1:na
      fxa=c4p3*(xa(ia)^alf3)*exp(-b*xa(ia)^gma);
      fxb=c4p3*(xb(ia)^alf3)*exp(-b*xb(ia)^gma);
      fv(ia)=(fxa+fxb)*da/2;
      vol=pi6*dia(ia)^3;
      rho(ia)=fv(ia)/vol;
      xi0=xi0+pi6*rho(ia);
      xi1=xi1+pi6*rho(ia)*dia(ia);
      xi2=xi2+pi6*rho(ia)*dia(ia)^2;
      xi3=xi3+fv(ia);
    end

%   build matrix Ct (k=0)
    Ct=zeros(na,na);
    dn=1.0-xi3;
    for ia=1:na
      ni=rho(ia);
      Ri=dia(ia);
      for ja=1:na
        nj=rho(ja);
        Rj=dia(ja);
        tm1=(Ri+Rj)^3/dn;
        tm2=(3*Ri*Rj*(Ri^2+Rj^2)*xi2+3*(Ri*Rj)^2*(Ri+Rj)*xi1+9*(Ri*Rj)^2*xi2+(Ri*Rj)^3*xi0)/(dn^2);
        tm3=(9*(Ri*Rj*xi2)^2*(Ri+Rj)+6*(Ri*Rj)^3*xi1*xi2)/(dn^3);
        tm4=(9*Ri^3*Rj^3*xi2^3)/(dn^4);
        Ct(ia,ja)=-pi6*(tm1+tm2+tm3+tm4)*sqrt(ni*nj);
      end
    end

%   build matrix Ht (k=0)
    Ht=zeros(na,na);
    Ht=inv(I-Ct)-I;

%   build matrix H (k=0)
    H=zeros(na,na);
    for ia=1:na
      ni=rho(ia);
      for ja=1:na
        nj=rho(ja);
        H(ia,ja)=Ht(ia,ja)/(sqrt(ni*nj));
      end
    end

    tm2=zeros(na,1);
    for ia=1:na;
      for ja=1:na
        tm2(ia)=tm2(ia)+dia(ja)^3*rho(ja)*H(ja,ia)/8;
      end
    end

%   initial solution
    p2=ones(1,3);
    p2(3)=mu*(ftot-1);
    p2(2)=mu-4*mu*ftot-1;
    root2=sqrt(roots(p2));
    kef=k*root2(1);

%   iterative solution
    for iter=1:niter
      y=mu/((kef/k)^2+mu);
      D=1-ftot*y;
      tm1=0.0;
      for ia=1:na
        tm1=tm1+fv(ia)*(1+i*2*kef^3*y*(dia(ia)^3/8+tm2(ia))/3/D);
      end
      kef=sqrt(k^2+3*kef^2*y*tm1/D);
    end
    ner=2*imag(kef)/k;
    edc=real(kef^2/k^2);
  end
  fprintf(fp,'%9.6f %14.9f %14.9f \n',ftot,ner,edc);
end
fclose(fp);

% plot attenuation rate and effective permittivity
clear;
load qcacp.dat -ascii;
figure;
plot(qcacp(:,1),qcacp(:,2),'r:');
axis([0.0,0.4,0.0,8e-4]);
xlabel('Volume Fraction ');
ylabel('Attenuation Rate ');
legend('QCA-CP-PY');
figure;
plot(qcacp(:,1),qcacp(:,3),'r:');
axis([0.0,0.4,1.0,1.7]);
xlabel('Volume Fraction ');
ylabel('Effective Permittivity ');
legend('QCA-CP-PY');
