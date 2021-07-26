function [M,C] = weights_colormap(S,nin)
  % Set up color map and color axis for a given scalar weight function
  %
  % [M,C] = weights_colormap(S,nin)
  %
  % Inputs:
  %   S  #V by 1 list of vertex scalar values
  %   nin  number of isointervals
  % Outputs:
  %   M  nin by 3 list of colors
  %   C  color axis min and max
  %
  % Example:
  %   % S contains some scalar weights on (V,F) mesh
  %   [M,C] = weights_colormap(S,10);
  %   trisurf(F,V(:,1),V(:,2),S);
  %   colormap(M);
  %   caxis(C);
  %   

  if max(S) > 1 && min(S) <0
    M = [...
      [linspace(1,0,nin); linspace(0,0,nin);linspace(0.5,0.5,nin);]'; ...
      jet(nin);
      [max(linspace(0.5,-0.5,nin),0.25); linspace(0,0.1,nin);min(linspace(0,0.2,nin),0.1);]'; ...
      ];
    C = [-1 2];
  elseif min(S) < 0
    M = [...
      [linspace(1,0,nin); linspace(0,0,nin);linspace(0.5,0.5,nin);]'; ...
      jet(nin)];
    C = [-1 1];
  else
    M = jet(nin);
    C = [0 1];
  end
end
