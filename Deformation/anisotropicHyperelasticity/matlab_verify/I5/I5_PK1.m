function final = I5_PK1(U, Sigma, V, a)
  F = U * Sigma * V';
  final = 2 * F * (a * a');
end
