function [H] = Corotational_Hessian_sym(F)
  eigenvector_gallery_sym;

  % get the singular values in an order that is consistent
  % with the numbering in the paper
  sigma0 = F(1,1);
  sigma1 = F(2,2);
  
  IS = trace(F);
  kTerm = (IS - 2);
  
  H = 2 * sym(eye(4,4));
  H = H + 2 * ((kTerm - 2) / (sigma0 + sigma1)) * (t * t');
  H = H + 2 * (r * r');
end
