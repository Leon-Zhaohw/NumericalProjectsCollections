invSqrt2 = sym(1) / sqrt(sym(2));
invSqrt3 = sym(1) / sqrt(sym(3));

T0 = sym([0 -1 0; 1 0 0; 0 0 0]);
T0 = invSqrt2 * T0;
  
T1 = sym([0 0 0; 0 0 1; 0 -1 0]);
T1 = invSqrt2 * T1;

T2 = sym([0 0 1; 0 0 0; -1 0 0]);
T2 = invSqrt2 * T2;

L0 = sym([0 1 0; 1 0 0; 0 0 0]);
L0 = invSqrt2 * L0;
  
L1 = sym([0 0 0; 0 0 1; 0 1 0]);
L1 = invSqrt2 * L1;

L2 = sym([0 0 1; 0 0 0; 1 0 0]);
L2 = invSqrt2 * L2;

S0 = sym([1 0 0; 0 0 0; 0 0 0]);
S1 = sym([0 0 0; 0 1 0; 0 0 0]);
S2 = sym([0 0 0; 0 0 0; 0 0 1]);

R = invSqrt3 * sym([1 0 0; 0 1 0; 0 0 1]);

% get flattened versions of everybody
t0 = vec(T0);
t1 = vec(T1);
t2 = vec(T2);

l0 = vec(L0);
l1 = vec(L1);
l2 = vec(L2);

s0 = vec(S0);
s1 = vec(S1);
s2 = vec(S2);

r = vec(R);
