invSqrt2 = 1 / sqrt(2);
invSqrt3 = 1 / sqrt(3);

T = [0 -1; 1 0];
T = invSqrt2 * T;
  
L = [0 1; 1 0];
L = invSqrt2 * L;
  
P = [1 0; 0 -1];
P = invSqrt2 * P;

S0 = [1 0; 0 0];
S1 = [0 0; 0 1];

R = invSqrt2 * [1 0; 0 1];

% get flattened versions of everybody
t = vec(T);
l = vec(L);
p = vec(P);

s0 = vec(S0);
s1 = vec(S1);

r = vec(R);
