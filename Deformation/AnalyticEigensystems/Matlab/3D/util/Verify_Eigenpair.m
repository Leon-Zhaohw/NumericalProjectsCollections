function [] = Verify_Eigenpair(H, q, lambda)
  trial = simplify(H * q);
  diff = simplify(trial - lambda * q);
  diffNorm = norm(diff);

  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  isMatlab = ~isOctave;

  if (isMatlab)
    syms sigma0 sigma1 sigma2 real;
    syms I1 I2 I3 real;
    shorter = lambda;
    shorter = subs(shorter, sigma0 + sigma1 + sigma2, I1);
    shorter = subs(shorter, sigma0^2 + sigma1^2 + sigma2^2, I2);
    shorter = subs(shorter, sigma0 * sigma1 * sigma2, I3);
    shorter = simplify(shorter);
  end

  allZeros = isequal(diffNorm, sym(0));
  if (allZeros)
      if (isMatlab)
          fprintf('Symbolic eigenvalue expression: %s VERIFIED\n', shorter);
      else
          fprintf('Symbolic eigenvalue: VERIFIED\n');
      end
  else
      if (isMatlab)
          fprintf('Symbolic eigenvalue expression: %s FAILED\n', shorter);
      else
          fprintf('Symbolic eigenvalue: FAILED\n');
          trial
          lq = simplify(lambda * q);
          simplify(diff)
      end
  end
end
