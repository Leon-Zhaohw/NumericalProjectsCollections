function t = Verify_Anisotropic_PK1(Psi,PK1)
  % Pick your sigmas
  Sigma = [
   2 0 0,
   0 5 0,
   0 0 11];
  a = [1 2 3]';
  a = a / norm(a);

  % Pick some random rotations and axes
  angleU = 95.1413;
  U = axisAngle([1 0.5 0.25], angleU);
  
  angleV = 212.1818;
  V = axisAngle([0.5 2.125 3.25], angleV);
  
  % compose your F
  F = U * Sigma * V';

  PK1direct = PK1(U, Sigma, V, a);
  PK1numerical = zeros(3,3);
  eps = 1e-2;

  % track whether or not the differences are converging
  previousdiff = inf;
  failed = false;
  
  for e = 1:5
    F0 = F;
    P0 = Psi(U, Sigma, V, a);
    
    for y = 1:3
      for x = 1:3
        F1 = F;
        F1(x,y) = F1(x,y) + eps;
        [U1 Sigma1 V1] = svd_rv(F1);
        
        P1 = Psi(U1, Sigma1, V1, a);
        diff = (P1 - P0) / eps;
        
        PK1numerical(x,y) = diff;
      end
    end
  
    PK1diff = PK1direct - PK1numerical;
    diff = norm(PK1diff) / 9;
    
    % see if we stopped converging
    if (diff >= (previousdiff / 2))
        failed = true;
    end
    previousdiff = diff;
    fprintf('eps: %.10f \t diff: %.10f\n', eps, diff)
    eps = eps * 0.1;
  end
  if (failed)
      fprintf('PK1 convergence test ***FAILED***\n')
  else
      fprintf('PK1 convergence test ***PASSED***\n')
  end
end
