function [final] = gradientJ(F)
  f0 = F(:,1);
  f1 = F(:,2);
  f2 = F(:,3);
  final = [cross(f1,f2) cross(f2,f0) cross(f0,f1)];
end
