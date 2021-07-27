% Boilerplate to load up the symbolic package in Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave)
    pkg load symbolic;
end

% Declare the invariants and constants
syms I5 mu real;

% Build the Anisotropic StVK energy
Psi_AS = (mu / 2) * (I5 - 1)^2;
[lambdas_AS] = Get_Analytic_Anisotropic_Eigenvalues(Psi_AS)

% Build the Anisotropic Sqrt energy
Psi_Sqrt = (mu / 2) * (sqrt(I5) - 1)^2;
[lambdas_Sqrt] = Get_Analytic_Anisotropic_Eigenvalues(Psi_Sqrt)

