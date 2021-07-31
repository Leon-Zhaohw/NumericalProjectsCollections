 function [varargout]=legslb(n)

%  x=legslb(n) returns n Legendre-Gauss-Lobatto points with x(1)=-1, x(n)=1
%  [x,w]= legslb(n) returns n Legendre-Gauss-Lobatto points and weights
%  Newton iteration method is used for computing nodes
%  See Page 99 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: lepoly() 
%  Last modified on August 30, 2011

 % Compute the initial guess of the interior LGL points
  nn=n-1; thetak=(4*[1:nn]-1)*pi/(4*nn+2);
  sigmak=-(1-(nn-1)/(8*nn^3)-(39-28./sin(thetak).^2)/(384*nn^4)).*cos(thetak);
  ze=(sigmak(1:nn-1)+sigmak(2:nn))/2;  
  ep=eps*10;                            % error tolerance for stopping iteration
  ze1=ze+ep+1;
 
 while max(abs(ze1-ze))>=ep,            % Newton's iteration procedure
      ze1=ze;
      [dy,y]=lepoly(nn,ze);
      ze=ze-(1-ze.*ze).*dy./(2*ze.*dy-nn*(nn+1)*y);  % see Page 99 of the book
 end;                                   % around 6 iterations are required for n=100
    varargout{1}=[-1,ze,1]';
 if nargout==1, return; end;
   
 % Use the weight expression (3.188) to compute the weights
  varargout{2}=[2/(nn*(nn+1)),2./(nn*(nn+1)*y.^2),2/(nn*(nn+1))]'; 




