% Pick your sigmas
Sigma = [2 0,
         0 5];

% Pick some rotations
angleU = 0;
U = [cos(angleU) -sin(angleU);
     sin(angleU)  cos(angleU)];

angleV = 0;
V = [cos(angleV) -sin(angleV);
     sin(angleV)  cos(angleV)];

% compose your F
F = U * Sigma * V';

% compose your polar decomposition
R = U * V';
S = V * Sigma * V';

% get the rotation gradient
H = DRDF_Horrible(R,S)

% what's the eigendecomposition look like?
[Q Lambda] = eig(H)

% is the eigenvalue a rational number?
rats(Lambda(4,4))
