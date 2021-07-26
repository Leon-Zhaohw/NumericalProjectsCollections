function [] = Verify_Eigenpair(H, q, lambda)
  trial = simplify(H * q);
  diff = simplify(trial - lambda * q);
  diffNorm = norm(diff);

  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  isMatlab = ~isOctave;

  allZeros = isequal(diffNorm, sym(0));

  if (allZeros)
      if (isMatlab)
          fprintf('Symbolic eigenvalue expression: %s VERIFIED\n', lambda);
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
