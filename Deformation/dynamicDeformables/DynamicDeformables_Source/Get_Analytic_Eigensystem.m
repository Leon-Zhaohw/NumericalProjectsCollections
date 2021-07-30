function [Q lambdas] = Get_Analytic_Eigensystem(Psi)
  syms I1 I2 I3 real;
  syms sigma0 sigma1 sigma2 sigmai sigmaj sigmak real;
  lambdas = sym(zeros(9,1));

  % x-twist, y-twist, and z-twist
  lambdas(4) = (2 / (sigma1 + sigma2)) * diff(Psi, I1) + ...
                2 * diff(Psi, I2) + sigma0 * diff(Psi, I3);
  lambdas(6) = (2 / (sigma0 + sigma2)) * diff(Psi, I1) + ...
                2 * diff(Psi, I2) + sigma1 * diff(Psi, I3);
  lambdas(5) = (2 / (sigma0 + sigma1)) * diff(Psi, I1) + ...
                2 * diff(Psi, I2) + sigma2 * diff(Psi, I3);

  % x-flip, y-flip, and z-flip
  lambdas(7) = 2 * diff(Psi, I2) - sigma2 * diff(Psi, I3);
  lambdas(8) = 2 * diff(Psi, I2) - sigma0 * diff(Psi, I3);
  lambdas(9) = 2 * diff(Psi, I2) - sigma1 * diff(Psi, I3);

  % x-scale, y-scale and z-scale
  lambdaScale = 2 * diff(Psi, I2) + diff(Psi, I1, 2) + ...
                4 * sigmai^2 * diff(Psi, I2, 2) + ...
                sigmaj^2 * sigmak^2 * diff(Psi, I3, 2) + ...
                4 * sigmai * diff(diff(Psi, I1), I2) + ...
                4 * I3 * diff(diff(Psi, I2), I3) + ...
                2 * sigmaj * sigmak * diff(diff(Psi, I3), I1);
  lambdas(1) = subs(lambdaScale, {sigmai, sigmaj, sigmak}, ... 
                                 {sigma0, sigma1, sigma2});
  lambdas(2) = subs(lambdaScale, {sigmai, sigmaj, sigmak}, ... 
                                 {sigma1, sigma0, sigma2});
  lambdas(3) = subs(lambdaScale, {sigmai, sigmaj, sigmak}, ... 
                                 {sigma2, sigma0, sigma1});

  % get compact expressions for the eigenvalues
  lambdas = subs(lambdas, I1, sigma0 + sigma1 + sigma2);
  lambdas = subs(lambdas, I2, sigma0^2 + sigma1^2 + sigma2^2);
  lambdas = subs(lambdas, I3, sigma0 * sigma1 * sigma2);
  lambdas = simplify(lambdas);
  
  % build the symbolic eigenvector gallery
  eigenvector_gallery_sym;

  % populate the eigenvectors
  Q = [s0 s1 s2 t0 t1 t2 l0 l1 l2];
end
