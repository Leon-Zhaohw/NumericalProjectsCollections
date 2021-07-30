function [final] = DRDF_Horrible(R,S)
  H = zeros(4,4);
  for index = 1:4
    i = mod(index - 1, 2);
    j = floor((index - 1) / 2);

    % promote it by a dimension  
    i = i + 1;
    j = j + 1;
    index3 = i + (j - 1) * 3;

    R3 = eye(3,3);
    S3 = eye(3,3);
    R3(1:2, 1:2) = R;
    S3(1:2, 1:2) = S;
  
    [DR3, DS3] = DRDF_Column(R3,S3,index3);
  
    DR = DR3(1:2,1:2);
    DS = DS3(1:2,1:2);
    
    column = DR;
    H(:,i) = reshape(column, 4, 1);
  end
  final = H;
endfunction
