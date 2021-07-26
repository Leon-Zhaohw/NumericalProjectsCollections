function [DF] = DFDF(index)
  i = mod(index - 1, 3);
  j = floor((index - 1) / 3);
  
  % Matlab indexing
  i = i + 1;
  j = j + 1;

  DF = zeros(3,3);
  DF(i,j) = 1;
end
