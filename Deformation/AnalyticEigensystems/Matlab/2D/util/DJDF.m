function [final] = DJDF(F)
  final = [F(2,2) -F(2,1);
           -F(1,2) F(1,1)];
end
