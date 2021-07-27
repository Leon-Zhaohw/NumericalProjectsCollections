      function[Ih] =  Interp_mat(xo,xi)
%
%     Compute the interpolation matrix from xi to xo
%
      
      no = length(xo);
      ni = length(xi);
      Ih = zeros(ni,no);
      w  = zeros(ni,2);
      for i=1:no;
         w = fd_weights_full(xo(i),xi,1);
         Ih(:,i) = w(:,1);
      end;
      Ih = Ih';
