function [H] = hessianJ(F)
  zero3 = zeros(3,3);
  f0 = F(:,1);
  f1 = F(:,2);
  f2 = F(:,3);
  f0hat = crossMatrix(f0); 
  f1hat = crossMatrix(f1); 
  f2hat = crossMatrix(f2); 
  H = [ zero3 -f2hat  f1hat;
        f2hat  zero3 -f0hat;
       -f1hat  f0hat  zero3];
end

function [final] = crossMatrix(f)
  final = [    0 -f(3)  f(2);
            f(3)     0 -f(1);
           -f(2)  f(1)     0];
end
