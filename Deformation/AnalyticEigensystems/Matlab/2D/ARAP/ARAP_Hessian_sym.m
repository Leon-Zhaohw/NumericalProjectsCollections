function [H] = ARAP_Hessian_sym(F)
  eigenvector_gallery_sym;
    
  % get the singular values in an order that is consistent
  % with the numbering in the paper
  s0 = F(1,1);
  s1 = F(2,2);
  
  H = 2 * sym(eye(4,4));
  H = H - (4 / (s0 + s1)) * (t * t');
end
