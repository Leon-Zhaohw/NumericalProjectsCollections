function [DR, DS] = DRDF_Column(R,S,index)
  i = mod(index - 1, 3);
  j = floor((index - 1) / 3);
  
  % Matlab indexing
  i = i + 1;
  j = j + 1;

  eiej = zeros(3,3);
  eiej(i,j) = 1;  

  G = (eye(3,3) * trace(S) - S) * R';
  Rij = R' * eiej;
  
  % extract the skew vector
  Rijsym = (Rij - Rij') * 0.5;
  skew = zeros(3,1);
  
  skew(1) = -Rijsym(2,3);
  skew(2) =  Rijsym(1,3);
  skew(3) = -Rijsym(1,2);
  skew = 2 * skew;
  
  omega = (G^-1) * skew;
  
  cross = [        0 -omega(3)  omega(2); 
            omega(3)         0 -omega(1);
           -omega(2)  omega(1)         0];
  DR = cross * R;
  DS = R' * (eiej - DR * S);
end
