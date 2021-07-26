function t = Verify_Hessian(PK1,Hessian)
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

  % track whether or not the differences are converging
  previousdiff = inf;
  failed = false;  
  
  Hdirect = Hessian(F);
  Hnumerical = zeros(4,4);
  eps = 1e-2;
  for e = 1:5
    F0 = F;
    P0 = PK1(F0);
    
    i = 1;
    for y = 1:2
      for x = 1:2
        F1 = F;
        F1(x,y) = F1(x,y) + eps;
        
        P1 = PK1(F1);
        diff = (P1 - P0) / eps;
        
        Hnumerical(:,i) = reshape(diff, 4,1);
        i = i + 1;
      end
    end
  
    Hdiff = Hdirect - Hnumerical;
    diff = norm(Hdiff) / 16;
    
    % see if we stopped converging
    if (diff >= (previousdiff / 2))
        failed = true;
    end
    previousdiff = diff;
    
    fprintf('eps: %.10f \t diff: %.10f\n', eps, diff);
    eps = eps * 0.1;
  end
  if (failed)
      fprintf('Hessian convergence test ***FAILED***\n')
  else
      fprintf('Hessian convergence test ***PASSED***\n')
  end  
end
