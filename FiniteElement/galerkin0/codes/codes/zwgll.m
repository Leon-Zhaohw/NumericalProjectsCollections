      function [z,w] = zwgll(p);

%
%
% computes the p+1 Gauss-Lobatto-Legendre nodes z on [-1,1]
% i.e. the zeros of the first derivative of the Legendre polynomial
% of degree p plus -1 and 1
%
% and the p+1 weights w
%


n = p+1;

z(1:n)=0;
w(1:n)=0;

z(1)=-1;
z(n)= 1;


if p>1,
if p==2 
   z(2)=0; 
else
  M=zeros(p-1,p-1);
  for i=1:p-2,
    M(i,i+1)=(1/2)*sqrt((i*(i+2))/((i+1/2)*(i+3/2)));
    M(i+1,i)=M(i,i+1);
  end;

  [V,D]=eig(M);
  z(2:p)=sort(eig(M));
end;
end;

%compute the weights w

w(1)=2/(p*(n));
w(n)=w(1);

for i=2:p,

  x=z(i);

  z0=1;
  z1=x;
  for j=1:p-1,
    z2=x.*z1*(2*j+1)/(j+1)-z0*j/(j+1);
    z0=z1;
    z1=z2;
  end;
  w(i)=2/(p*(n)*z2*z2);

end;

z=z';
w=w';
