function final = I5_Psi(U, Sigma, V, a)
  F = U * Sigma * V';
  final = a' * F' * F * a;
end
