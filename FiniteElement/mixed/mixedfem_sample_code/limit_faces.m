function [LF] = limit_faces(F,L,exclusive)
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, NYU
  %

  % limit faces F to L 
  if(exist('exclusive') && exclusive)
    % all indices must be in F
    LF = F(ismember(F(:,1),L) & ismember(F(:,2),L) & ismember(F(:,3),L),:);
  else
    % at least one index must be in F
    LF = F(ismember(F(:,1),L) | ismember(F(:,2),L) | ismember(F(:,3),L),:);
  end
end
