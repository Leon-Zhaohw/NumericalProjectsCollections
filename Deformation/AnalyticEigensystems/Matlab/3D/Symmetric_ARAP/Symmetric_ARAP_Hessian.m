function [H] = Symmetric_ARAP_Hessian(F)
  [U Sigma V] = svd(F);

  djdf = DJDF(F);
  g = vec(djdf);

  J = det(F);
  C = F' * F;
  IC = trace(C);
  IIC = trace(C * C);
  IIStarC = 0.5 * (IC^2 - IIC);

  dIIStarC = 2 * IC * F - 2 * F * F' * F;
  t = vec(dIIStarC);
  d2IIStarC = IIC_Star_Hessian(F);
  d2J = HessianJ(F);

  H = 2 * eye(9,9);
  H = H - (2 / J^3) * (g * t' + t * g');
  H = H + (6 * IIStarC / J^4) * (g * g');
  H = H + (1 / J^2) * d2IIStarC;
  H = H - (2 * IIStarC / J^3) * d2J;

  middle0 = [ 0 1 0;
             -1 0 0;
              0 0 0];
  middle1 = [ 0  0 0;
              0  0 1;
              0 -1 0];
  middle2 = [ 0 0 1;
              0 0 0;
             -1 0 0];

  Q0 = (1 / sqrt(2)) * U * middle0 * V';
  Q1 = (1 / sqrt(2)) * U * middle1 * V';
  Q2 = (1 / sqrt(2)) * U * middle2 * V';

  q0 = vec(Q0);
  q1 = vec(Q1);
  q2 = vec(Q2);

  s0 = Sigma(1,1);
  s1 = Sigma(2,2);
  s2 = Sigma(3,3);

  DRDF =        (2 / (s0 + s1)) * (q0 * q0');
  DRDF = DRDF + (2 / (s1 + s2)) * (q1 * q1');
  DRDF = DRDF + (2 / (s0 + s2)) * (q2 * q2');

  R = U * V';
  S = V * Sigma * V';
  r = vec(R);
  f = vec(F);

  IS = trace(S);
  IIS = trace(S' * S);
  IIStarS = 0.5 * (IS^2 - IIS);

  newH = (1 + IS / J) * DRDF + (1 / J) * r * r';
  newH = newH - (IS / J^2) * (g * r' + r * g');
  newH = newH + (2.0 * IIStarS / J^3) * (g * g');
  newH = newH - (IIStarS / J^2) * d2J;
  newH = newH + (1 / J^2) * (g * f' + f * g'); 
  newH = newH - (1/J) * eye(9,9);
  H = H - 2.0 * newH;
  
  H = H / 2;
end
