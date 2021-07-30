% Boilerplate to load up the symbolic package in Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
    pkg load symbolic;
end

% Declare the invariants and constants
syms I1 I2 I3 mu lambda real;

% Build the SNH energy
% Psi = (mu / 2) * (I2 - 3) - mu * (I3 - 1) + (lambda / 2) * (I3 - 1)^2;
% Build the ARAP energy
Psi = I2 - 2 * I1 + 3

% Get the analytic twist and flip eigenvalues, and
% maybe the stretching values too, if you hti the jackpot
[lambdas] = Get_Analytic_Eigenvalues(Psi);

% Get the stretching matrix
A = Get_Stretching_System(Psi);

% Did we hit the Simple Eigenvalue Jackpot?
jackpot = A(1,2) + A(1,3) + A(2,3);
if (jackpot == sym(0))
  printf('You HIT the Jackpot!\n');
  printf('No stretching matrix is needed.\n');
  printf('Here are your eigenvalues:\n');
  lambdas
else
  printf('You MISSED the Jackpot!\n');
  printf('The stretching eigenvalues are the three\n');
  printf('eigenvalues from this matrix:\n');
  A
  printf('The last six are these:\n');
  lambdas(4:9)
end
