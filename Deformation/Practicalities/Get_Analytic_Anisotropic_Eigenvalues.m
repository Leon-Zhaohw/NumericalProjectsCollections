function [lambdas] = Get_Analytic_Eigenvalues(Psi)
  syms mu I5 real;
  lambdas = sym(zeros(3,1));

  firstI5 = diff(Psi, I5);
  secondI5 = diff(firstI5, I5);

  lambdas(1) = 2 * (firstI5 + 2 * I5 * secondI5);
  lambdas(2) = 2 * firstI5;
  lambdas(3) = 2 * firstI5;

  lambdas = simplify(lambdas);
end
