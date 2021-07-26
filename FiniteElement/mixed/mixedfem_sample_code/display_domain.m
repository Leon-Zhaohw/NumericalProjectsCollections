function display_domain(F,V,Omega,N0,N1,N2,outside_region_of_interest,phong)
%
% See corresponding paper: "Mixed finite elements for variational surface
% modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
% 2010
%
% Copyright 2010, Alec Jacobson, NYU
%

  trisurf( ...
    limit_faces(F,Omega), ...
    V(:,1),V(:,2),V(:,3), ...
    'FaceColor',[1.0,0.9,0.1]);
  hold on;
  trisurf( ...
    limit_faces(F,[N0, N1, N2, outside_region_of_interest],1), ...
    V(:,1),V(:,2),V(:,3), ...
    'FaceColor',[0.2,0.1,0.6]);
  hold on;
  % Display rows in the exterior
  plot3(V(N0,1),V(N0,2),V(N0,3), ...
    'o','MarkerEdgeColor',[0,0,0], 'MarkerFaceColor',[0.5,0.1,0.1]);
  hold on;
  plot3(V(N1,1),V(N1,2),V(N1,3), ...
    's','MarkerEdgeColor',[0,0,0], 'MarkerFaceColor',[0.1,0.1,0.5]);
  hold on;
  plot3(V(N2,1),V(N2,2),V(N2,3), ...
    'd','MarkerEdgeColor',[0,0,0], 'MarkerFaceColor',[0.1,0.5,0.1]);
  view(2);
  legend('Omega','Exterior','N0','N1','N2');
  
  if(exist('phong') & phong)
    light_positions = {
      [-1.0,-1.0,2.0], ...
      [0.0,-1.0,0.0], ...
      [0.0,1.0,0.0], ...
      };
    %phong lighting
    for l=1:size(light_positions,2)
      light('Position',light_positions{l},'Style','infinite');
    end
    lighting phong;
  end
end
