function [final] = DJDF(F)
  f0 = F(:,1);
  f1 = F(:,2);
  f2 = F(:,3);
  final = zeros(3,3);
  final(:,1) = cross(f1,f2);
  final(:,2) = cross(f2,f0);
  final(:,3) = cross(f0,f1);
end
