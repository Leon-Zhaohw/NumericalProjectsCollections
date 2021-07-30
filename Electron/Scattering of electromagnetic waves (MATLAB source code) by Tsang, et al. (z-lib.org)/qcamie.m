function K_eff=qcamie(freq,epsilon_p,f,k0a,n_max)
%QCAMIE computes the effective propagation constant using the quasi-crystalline %approximation (QCA) for a medium consisting of densely distributed Mie
%scatterers.
%
%    K_eff=qcamie(freq,epsilon_p,f,nmax,ka)
%
%    INPUT:
%
%    freq=frequency in GHz	
%    epsilon_p=particle permittivity relative to homogeneous background
%    f=fractional volume of particles
%    k0a=size parameter (this can be a vector)
%    n_max=maximum spherical multipole used
%
%    OUTPUT:
%
%    K_eff=complex number (per cm) which denote the effective propagation 
%          constant at each k0a
%
% -- Part of the Electromagnetic Wave MATLAB Library (EWML)--
%    <http://www.emwave.com/>

% Original: by Chite Chen, November 1998

tol=1e-14;                         % error tolerance of det(T)
lambda=30/freq;                    % wavelength in cm
k=2*pi/lambda;	                 % wavenumber in 1/cm
na=max(size(k0a));

% Read pair function for given f for the integration limit of Mp
load pair.dat;
r_b=pair(:,1);
gg=pair(:,2);

for ia=1:na,
  ka=k0a(ia);
  a=ka/k;
  b=2*a;
  kpa=ka*sqrt(epsilon_p);
  no=6*f/(pi*b^3);

% First initial guess is the solution for media with sparse concentration
% the exciting field approximately the same as the incident field
  FF=0;
  for nn=1:n_max,
    FF=FF+(2*nn+1)*(Tn_M(nn,ka,kpa)+Tn_N(nn,ka,kpa));
  end
  K_F=k-i*pi*no/k^2*FF;

% Second initial guess is the low frequency limit solution
% For Percus-Yevick pair function
  y=(epsilon_p-1)/(epsilon_p+2);
  K_low=sqrt(k^2+3*f*k^2*y/(1-f*y)*(1+i*2/3*(ka)^3*y*(1-f)^4/((1-f*y)*(1+2*f)^2)));

% Third initial guess
  K_low_real=real(K_low);

  r=r_b*b;                          % resize r;

  x1=K_F;
  x2=real(K_low);
  x3=K_low;
  T1=SysEqu(n_max,k,x1,ka,kpa,b,no,r,gg);
  T2=SysEqu(n_max,k,x2,ka,kpa,b,no,r,gg);
  T3=SysEqu(n_max,k,x3,ka,kpa,b,no,r,gg);
  f1=det(T1);
  f2=det(T2);
  f3=det(T3);

% Rearrange the order so that abs(f(i))<= abs(f(i-1)) <= abs(f(i-2))
  [x1 x2 x3 f1 f2 f3]=Rearrange(x1,x2,x3,f1,f2,f3);
  if abs(f3)<tol,
    K_eff(ia)=x3;
  else
% Using Muller's Approach to calculate the new guesses
    n_iter=0;                       % number of iteration
    while (abs(f3)>tol)&(n_iter<=30),
      x_new=Muller(x1,x2,x3,f1,f2,f3);
      x1=x2;
      x2=x3;
      x3=x_new;
      T1=T2;
      T2=T3;
      T3=SysEqu(n_max,k,x3,ka,kpa,b,no,r,gg);
      f1=f2;
      f2=f3;
      f3=det(T3);
      [x1 x2 x3 f1 f2 f3]=Rearrange(x1,x2,x3,f1,f2,f3);
      n_iter=n_iter+1;
    end
    K_eff(ia)=x3;
  end
  K_eff
end

figure(1)
plot(k0a,k./real(K_eff));
grid on;
axis([0 2.5 0.8 1])
xlabel('ka')
ylabel('Phase Velocity')
title('Normalized phase velocity k/K_r as function of ka')
%
figure(2)
semilogy(k0a,2*imag(K_eff)./real(K_eff));
grid on;
axis([0 2.5 1e-4 1])
xlabel('ka')
ylabel('Loss Tangent')
title('Effective loss tangent 2K_i/K_r as function of ka')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Tn_M.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Tn_M(n,ka,kpa)
% compute T-matrix elements for vector spherical waves M_mn

numerator1=sbesselj(n,kpa).*(sbesselj(n,ka)+ka*sbesselj_p(n,ka));
numerator2=sbesselj(n,ka).*(sbesselj(n,kpa)+kpa*sbesselj_p(n,kpa));
denominator1=sbesselj(n,kpa).*(sbesselh(n,ka)+ka*sbesselh_p(n,ka));
denominator2=sbesselh(n,ka).*(sbesselj(n,kpa)+kpa*sbesselj_p(n,kpa));
output=-(numerator1-numerator2)./(denominator1-denominator2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Tn_N.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Tn_N(n,ka,kpa)
% compute T-matrix elements for vector spherical waves N_mn

numerator1=(kpa)^2*sbesselj(n,kpa).*(sbesselj(n,ka)+ka*sbesselj_p(n,ka));
numerator2=(ka)^2*sbesselj(n,ka).*(sbesselj(n,kpa)+kpa*sbesselj_p(n,kpa));
denominator1=(kpa)^2*sbesselj(n,kpa).*(sbesselh(n,ka)+ka*sbesselh_p(n,ka));
denominator2=(ka)^2*sbesselh(n,ka).*(sbesselj(n,kpa)+kpa*sbesselj_p(n,kpa));
output=-(numerator1-numerator2)./(denominator1-denominator2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Clebsch.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Clebsch(j1,j2,j3,m1,m2,m)
% calculate Clebsch-Gordan coefficients using the formula in
% Abramowitz and Stegun

FF=0;
if (j1 < abs(m1)) | (j2 < abs(m2)) | (j3 < abs(m)),
  output=0;
% sprintf('Condition (j1>=|m1|, j2 >=|m2| and j >=|m|) does NOT obey')
  break
elseif ((j3 > j1+j2) | j3 < abs(j1-j2)),
  output=0;
% sprintf('Condition ( |j1-j2|<=j<=j1+j2 ) does NOT obey')
  break
end

if (m1+m2)~= m,
  output=0;
else
  term1=1/2*(faclog(j1+j2-j3)+faclog(j3+j1-j2)+faclog(j3+j2-j1) ...
       +log(2*j3+1)-faclog(j3+j1+j2+1)+faclog(j1+m1)+faclog(j1-m1) ...
       +faclog(j2+m2)+faclog(j2-m2)+faclog(j3+m)+faclog(j3-m));
% term1 includes the terms which are not related to k
% determine the range for k - the factorial cannot be negative
  upperlimit=min([j1+j2-j3 j1-m1 j2+m2]);
  lowerlimit=abs(min([j3-j2+m1 j3-j1-m2 0]));
  for k=lowerlimit:upperlimit,
    term2=-(faclog(k)+faclog(j1+j2-j3-k)+faclog(j1-m1-k) ...
         +faclog(j2+m2-k)+faclog(j3-j2+m1+k)+faclog(j3-j1-m2+k));
    FF=FF+(-1)^k*exp(term2);
  end
  output=FF*exp(term1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             wigner.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=wigner(j1,j2,j3,m1,m2,m)
% calculate Wigner 3-j symbol

output=(-1)^(j1-j2-m)*(2*j3+1)^(-1/2)*Clebsch(j1,j2,j3,m1,m2,-m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            factorial.m                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=factorial(n)
% compute the factorial of n

product=1;
if n==0,
  output=1;
elseif (n < 0),
  sprintf('n cannot be negative')
  break
else
  for n_index=1:n,
    product=n_index*product;
  end
  output=product;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             faclog.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=faclog(n)
% natural logarithm of factorial n (log(n!)) preventing overflow

nn=0;
if (n==0),
  output=0;                         % log(1) = 0;
elseif (n<0),
  sprintf('n cannot be negative')
  break
else
  for n_index=1:n,
     nn=nn+log(n_index);
  end
  output=nn;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            sbesselj.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=sbesselj(n,arg)
% spherical Bessel function of order n
% allow arguement to be an array, but only allow a single order n

output=sqrt(pi./(2*arg)).*besselj(n+1/2,arg);

% special handle for arguement = 0
[xx]=find (arg==0);
if (n==0),
  output([xx])=1;
else
  output([xx])=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           sbesselj_p.m                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=sbesselj_p(n,arg)
% derivative of spherical Bessel function of order n

output=1/(2*n+1)*(n*sbesselj(n-1,arg)-(n+1)*sbesselj(n+1,arg));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            sbesselh.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=sbesselh(n,arg)
% spherical Hankel function of the first kind of order n
% allow arguement to be an array
% singular at arguement = 0

output = sqrt(pi./(2*arg)).*besselh(n+1/2,1,arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           sbesselh_p.m                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=sbesselh_p(n,arg)
% derivative of spherical Hankel function of order n

output=1/(2*n+1)*(n*sbesselh(n-1,arg)-(n+1)*sbesselh(n+1,arg));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             a_mnuvp.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=a_mnuvp(m,n,u,v,p)
% coefficient a(

output=(-1)^(m+u)*(2*p+1)*sqrt((factorial(n+m)*factorial(v+u)*factorial(p-m-u)) ...
      /(factorial(n-m)*factorial(v-u)*factorial(p+m+u)))*wigner(n,v,p,m,u,-(m+u)) ...
      *wigner(n,v,p,0,0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            a_mnuvpq.m                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=a_mnuvpq(m,n,u,v,p,q)
% coefficient a(

output=(-1)^(m+u)*(2*p+1)*sqrt((factorial(n+m)*factorial(v+u)*factorial(p-m-u)) ...
      /(factorial(n-m)*factorial(v-u)*factorial(p+m+u)))*wigner(n,v,p,m,u,-(m+u)) ...
      *wigner(n,v,q,0,0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              A_nvp.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=A_nvp(n,v,p)
% coefficient A(

output=1/(n*(n+1)*(2*v+1))*(2*v*(v+1)*(2*v+1)+(v+1)*(n+v-p)*(n+p-v+1) ...
      -v*(n+v+p+2)*(v+p-n+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              B_nvp.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=B_nvp(n,v,p)
% coefficient B(

output=1/(n*(n+1))*sqrt((n+v+p+1)*(v+p-n)*(n+p-v)*(n+v-p+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Lp.m                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Lp(p,k,Keff,b)
% coefficient L

output=-b^2/(Keff^2-k^2)*(k*sbesselh_p(p,k*b)*sbesselj(p,Keff*b) ...
      -Keff*sbesselh(p,k*b)*sbesselj_p(p,Keff*b));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Mp2.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Mp2(p,k,Keff,r,gg,b)
% Pair function is independent of p
% exclude those for r less than b
% find the cutoff point

nr=max(size(r));
r_cutoff=min(find(r>=b));
r=r(r_cutoff:nr);
gg=gg(r_cutoff:nr);
h=gg-1;
hp=sbesselh(p,k*r);
jp=sbesselj(p,Keff*r);
yy=r.^2.*h.*hp.*jp;
output=trapz(r,yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Muller.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=Muller(x1,x2,x3,f1,f2,f3)
% Muller's approach

lambda_i=(x3-x2)/(x2-x1);
delta_i=1+lambda_i;
c_i=f1*lambda_i^2-f2*delta_i^2+f3*(lambda_i+delta_i);
den1=c_i+sqrt(c_i^2-4*f3*delta_i*lambda_i*(f1*lambda_i-f2*delta_i+f3));
den2=c_i-sqrt(c_i^2-4*f3*delta_i*lambda_i*(f1*lambda_i-f2*delta_i+f3));
if abs(den1)>=abs(den2),
  denominator=den1;
else
  denominator=den2;
end
output=x3+(x3-x2)*(-2*f3*delta_i)/denominator;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             SysEqu.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=SysEqu(n_max,k,Keff,ka,kpa,b,no,r,gg)
% calculate Mp+Lp and store as an array.
% abs(n-v) <=p <=(n+v)
% the upper bound and lower bound are 0 and 2*n_max

for v=1:n_max,
  for n=1:n_max,
    FF1=0;
    FF2=0;
    for p=abs(n-v):(n+v);
      MLSUM=Lp(p,k,Keff,b)+Mp2(p,k,Keff,r,gg,b);
      aA=a_mnuvp(1,n,-1,v,p)*A_nvp(n,v,p);
      aB=a_mnuvpq(1,n,-1,v,p,p-1)*B_nvp(n,v,p);
      FF1=FF1+MLSUM*aA;
      FF2=FF2+MLSUM*aB;
    end
    MM(v,n)=-2*pi*no*(2*n+1)*Tn_M(n,ka,kpa)*FF1;
    MN(v,n)=-2*pi*no*(2*n+1)*Tn_N(n,ka,kpa)*FF2;
    NM(v,n)=-2*pi*no*(2*n+1)*Tn_M(n,ka,kpa)*FF2;
    NN(v,n)=-2*pi*no*(2*n+1)*Tn_N(n,ka,kpa)*FF1;
  end
end
TT=[MM MN; NM NN];
output=eye(2*n_max)-TT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Rearrange.m                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,x2,x3,f1,f2,f3]=Rearrange(x1,x2,x3,f1,f2,f3)
% random input order, rearrange the order so that f1 >= f2 >= f3

if (abs(f1)<abs(f2)),
  tmp=f2;
  tmpX=x2;
  f2=f1;
  x2=x1;
  f1=tmp;
  x1=tmpX;              % if (f1<f2) interchange f1 and f2 so f1 >= f2
end

if (abs(f2)<abs(f3)),
  tmp=f3;
  tmpX=x3;
  f3=f2;
  x3=x2;
  f2=tmp;
  x2=tmpX;
end

if (abs(f1)<abs(f2)),
  tmp=f2;
  tmpX=x2;
  f2=f1;
  x2=x1;
  f1=tmp;
  x1=tmpX;
end
output=[x1 x2 x3 f1 f2 f3];