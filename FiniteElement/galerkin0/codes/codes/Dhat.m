      function[Dh] =  Dhat(x)
%
%                                                 ^
%     Compute the interpolatory derivative matrix D_ij associated
%     with nodes x_j such that
%
%                ^
%            w = D*u   
%            -     -
%
%     returns the derivative of u at the points x_i.
%
%
      
      n1 = length(x);
      w  = zeros(n1,2);
      Dh = zeros(n1,n1);
      for i=1:n1;
         w = fd_weights_full(x(i),x,1);
         Dh(:,i) = w(:,2);
      end;
      Dh = Dh';
