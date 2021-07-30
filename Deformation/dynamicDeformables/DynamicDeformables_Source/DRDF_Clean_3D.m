function [H] = DRDF_Clean_3D(U, Sigma, V)
  % get the twist modes
  T0 = [0 -1  0; 
        1  0  0; 
        0  0  0];
  T0 = (1 / sqrt(2)) * U * T0 * V';
  
  T1 = [0  0  0; 
        0  0  1; 
        0 -1  0];
  T1 = (1 / sqrt(2)) * U * T1 * V';

  T2 = [ 0  0  1; 
         0  0  0; 
        -1  0  0];
  T2 = (1 / sqrt(2)) * U * T2 * V';

  % get the flattened versions
  t0 = vec(T0);
  t1 = vec(T1);
  t2 = vec(T2);

  % get the singular values
  sx = Sigma(1,1);
  sy = Sigma(2,2);
  sz = Sigma(3,3);
  
  H = (2 / (sx + sy)) * (t0 * t0');
  H = H + (2 / (sy + sz)) * (t1 * t1');
  H = H + (2 / (sx + sz)) * (t2 * t2');
end
