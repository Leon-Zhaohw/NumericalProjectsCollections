function [out] = Simplify_Invariants(in)
  syms I1 I2 I3 sigma0 sigma1 sigma2 real;
  out = subs(in, I1, sigma0 + sigma1 + sigma2);
  out = subs(out, I2, sigma0^2 + sigma1^2 + sigma2^2);
  out = subs(out, I3, sigma0 * sigma1 * sigma2);
  out = simplify(out);
  out = subs(out, sigma0 + sigma1 + sigma2, I1);
  out = subs(out, sigma0^2 + sigma1^2 + sigma2^2, I2);
  out = subs(out, sigma0 * sigma1 * sigma2, I3);
  out = simplify(out);
end
