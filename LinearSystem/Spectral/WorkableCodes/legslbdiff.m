% D=legslbdiff(n,x) returns the first-order differentiation matrix of size
% n by n, associated with the Legendre-Gauss-Lobatto points x, which may be computed by 
% x=legslb(n) or x=legslbndm(n). Note: x(1)=-1 and x(n)=1.
% See Page 110 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: lepoly() 
% Last modified on August 31, 2011

function D=legslbdiff(n,x)
if n==0, D=[]; return; end;
xx=x;y=lepoly(n-1,xx); nx=size(x); 
if nx(2)>nx(1), y=y'; xx=x'; end;  %% y is a column vector of L_{n-1}(x_k)
  D=(xx./y)*y'-(1./y)*(xx.*y)';  %% compute L_{n-1}(x_j) (x_k-x_j)/L_{n-1}(x_k);     
                                 % 1/d_{kj} for k not= j (see (3.203)) 
  D=D+eye(n);                    % add the identity matrix so that 1./D can be operated                                     
  D=1./D; 
  D=D-eye(n); D(1,1)=-n*(n-1)/4; D(n,n)=-D(1,1);  %update the diagonal entries  
  return; 
 