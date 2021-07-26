function t = Verify_PK1(Psi,PK1)
  % Pick your sigmas
  Sigma = [
   2 0, 
   0 11];

  % Pick some random rotations and axes
  angleU = 95.1413;
  U = axisAngle2D(angleU);
  
  angleV = 212.1818;
  V = axisAngle2D(angleV);
  
  % compose your F
  F = U * Sigma * V';

  PK1direct = PK1(F);
  PK1numerical = zeros(2,2);
  eps = 1e-2;
  
  % track whether or not the differences are converging
  previousdiff = inf;
  failed = false;
  
  for e = 1:5
    F0 = F;
    P0 = Psi(F0);
    
    for y = 1:2
      for x = 1:2
        F1 = F;
        F1(x,y) = F1(x,y) + eps;
        
        P1 = Psi(F1);
        diff = (P1 - P0) / eps;
        
        PK1numerical(x,y) = diff;
      end
    end
  
    PK1diff = PK1direct - PK1numerical;
    diff = norm(PK1diff) / 4;
    
    % see if we stopped converging
    if (diff >= (previousdiff / 2))
        failed = true;
    end
    previousdiff = diff;
    
    fprintf('eps: %.10f \t diff: %.10f\n', eps, diff)
    eps = eps * 0.1;
  end
  %PK1direct
  %PK1numerical
  if (failed)
      fprintf('PK1 convergence test ***FAILED***\n')
  else
      fprintf('PK1 convergence test ***PASSED***\n')
  end
end
