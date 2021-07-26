isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
    pkg load symbolic;
end
addpath('./util')
addpath('./I5')

fprintf('============================================\n');
fprintf('Running numerical tests for I5\n');
fprintf('============================================\n');
fprintf('Verifying PK1 for I5\n')
Verify_Anisotropic_PK1(@I5_Psi, @I5_PK1)
fprintf('Verifying Hessian for I5\n')
Verify_Anisotropic_Hessian(@I5_PK1, @I5_Hessian)

fprintf('============================================\n');
fprintf('Running symbolic tests for I5\n');
fprintf('============================================\n');
U_sym = sym(eye(3,3));
V_sym = sym(eye(3,3));
syms sigma0 sigma1 sigma2 real;
sigma_sym = diag([sigma0 sigma1 sigma2]);
syms a0 a1 a2 real;
a_sym = [a0 a1 a2]';

% get the Hessian
H_sym = I5_Hessian_sym(U_sym, sigma_sym, V_sym, a_sym);

% verify the eigenvectors
zeros_sym = sym(zeros(3,3));
H_accum = sym(zeros(9,9));
for i = 1:3
  Q = zeros_sym;

  % try a different row
  Q(i,:) = a_sym';
  q = vec(Q);
  lambda = sym(2) * (a_sym' * a_sym);
  Verify_Eigenpair(H_sym, q, lambda);

  % accumulate the outer products
  H_accum = H_accum + simplify(lambda * (q * q'));
end

% verify that these three eigenvectors span the entire I5 Hessian
% i.e. that it is rank-three and these are the eigenvectors

% apply the assumption that the fiber direction is normalized
H_accum = subs(H_accum, (a0^2 + a1^2 + a2^2), sym(1));
diff = simplify(H_accum - H_sym);
diffNorm = norm(diff);
allZeros = isequal(diffNorm, sym(0));
if (allZeros)
  fprintf('Span of proposed eigensystem is complete: VERIFIED\n');
else
  fprintf('Span of proposed eigensystem is complete: FAILED\n');
  diff
end
