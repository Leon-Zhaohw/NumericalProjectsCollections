isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
    pkg load symbolic;
end

syms sigma0 sigma1 sigma2 real;
syms I1 I2 I3 mu lambda real;

%Psi = (mu / 8) * (8 * I1 * I3 + I2^2 + 2 * I1^2 * I2 - 4 * I2 - I1^4 + 6) + (lambda / 8) * (I2 - 3)^2;
Psi = I2 - 2 * I1 + 3
%Psi = (mu / 8) * (8 * I1 * I3 + I2^2 + 2 * I1^2 * I2 - 4 * I2 - I1^4 + 6)
%Psi = (lambda / 8) * (I2 - 3)^2;

A = Get_Stretching_System(Psi)

jackpot = A(1,2) + A(1,3) + A(2,3)

if (jackpot == sym(0))
  jackpotHit = 1
else
  jackpotHit = 0
end

[V] = Get_Analytic_Eigenvalues(Psi)

