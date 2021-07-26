invSqrt2 = sym(1) / sqrt(sym(2));
invSqrt3 = sym(1) / sqrt(sym(3));

T = sym([0 -1; 1 0]);
T = invSqrt2 * T;
  
L = sym([0 1; 1 0]);
L = invSqrt2 * L;
  
P = sym([1 0; 0 -1]);
P = invSqrt2 * P;

S0 = sym([1 0; 0 0]);
S1 = sym([0 0; 0 1]);

R = invSqrt2 * sym([1 0; 0 1]);

% get flattened versions of everybody
t = vec(T);
l = vec(L);
p = vec(P);

s0 = vec(S0);
s1 = vec(S1);

r = vec(R);
