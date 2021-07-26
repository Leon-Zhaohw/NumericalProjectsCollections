function [H] = Symmetric_ARAP_Hessian_sym(F)
  djdf = DJDF_sym(F);
  g = vec(djdf);

  J = det(F);
  C = F' * F;
  IC = trace(C);
  IIC = trace(C * C);
  IIStar = (IC^2 - IIC) / 2;

  dIIStar = 2 * IC * F - 2 * F * F' * F;
  t = vec(dIIStar);
  d2IIStar = IIC_Star_Hessian_sym(F);
  d2J = HessianJ(F);

  H = 2 * eye(9,9);
  H = H - (2 / J^3) * (g * t' + t * g');
  H = H + (6 * IIStar / J^4) * (g * g');
  H = H + (1 / J^2) * d2IIStar;
  H = H - (2 * IIStar / J^3) * d2J;

  middle0 = [ 0 1 0;
             -1 0 0;
              0 0 0];
  middle1 = [ 0  0 0;
              0  0 1;
              0 -1 0];
  middle2 = [ 0 0 1;
              0 0 0;
             -1 0 0];
  middle0 = sym(middle0);
  middle1 = sym(middle1);
  middle2 = sym(middle2);

  invSqrt = 1 / sqrt(sym(2));

  Q0 = middle0 * invSqrt;
  Q1 = middle1 * invSqrt;
  Q2 = middle2 * invSqrt;

  q0 = vec(Q0);
  q1 = vec(Q1);
  q2 = vec(Q2);

  s0 = F(1,1);
  s1 = F(2,2);
  s2 = F(3,3);

  DRDF =        (2 / (s0 + s1)) * (q0 * q0');
  DRDF = DRDF + (2 / (s1 + s2)) * (q1 * q1');
  DRDF = DRDF + (2 / (s0 + s2)) * (q2 * q2');

  R = sym(eye(3,3));
  r = vec(R);
  f = vec(F);

  S = F;
  IS = trace(S);
  IIS = trace(S' * S);
  IIStarS = (IS^2 - IIS) / 2;

  newH = (1 + IS / J) * DRDF + (1 / J) * r * r';
  newH = newH - (IS / J^2) * (g * r' + r * g');
  newH = newH + (2 * IIStarS / J^3) * (g * g');
  newH = newH - (IIStarS / J^2) * d2J;
  newH = newH + (1 / J^2) * (g * f' + f * g'); 
  newH = newH - (1/J) * sym(eye(9,9));
  H = H - 2 * newH;
  
  H = H / 2;
end
