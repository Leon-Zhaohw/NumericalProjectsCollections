function [] = Verify_Eigenpair(H, q, lambda)
  trial = simplify(H * q);
  diff = simplify(trial - lambda * q);
  diffNorm = norm(diff);

  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  isMatlab = ~isOctave;

  allZeros = isequal(diffNorm, sym(0));

  if (allZeros)
      if (isMatlab)
          syms sigma0 sigma1 real;
          syms I1 I2 I3 real;
          shorter = lambda;
          shorter = subs(shorter, sigma0 + sigma1, I1);
          shorter = subs(shorter, sigma0 * sigma1, I3);
          fprintf('Symbolic eigenvalue expression: %s VERIFIED\n', shorter);
      else
          fprintf('Symbolic eigenvalue: VERIFIED\n');
      end
  else
      if (isMatlab)
          fprintf('Symbolic eigenvalue expression: %s FAILED\n', lambda);
      else
          fprintf('Symbolic eigenvalue: FAILED\n');
      end
  end
end
