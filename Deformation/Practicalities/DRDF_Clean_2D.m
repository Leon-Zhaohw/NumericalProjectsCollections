function [final] = DRDF_Clean_2D(U, Sigma, V)
  Q0 = (1 / sqrt(2)) * U * [0 -1; 1 0] * V';
  q0 = vec(Q0);
  lambda0 = 2 / (Sigma(1,1) + Sigma(2,2));
  final = lambda0 * (q0 * q0');
endfunction
