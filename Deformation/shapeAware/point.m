function s = point(a,P,c,size_data)
  % POINT plot P with color c
  %
  % s = point(a,P,c,size_data)
  %
  % Inputs:
  %   a  axis to plot into
  %   P  #P by 3 list of points
  %   c  RGB color
  %   size_data size data of points
  % Outputs:
  %   s  handles to plotted points
  %

  offset = 1e-6;
  s1 = scatter3(a,P(:,1),P(:,2),offset + P(:,3),'o', ...
    'MarkerFaceColor','none','MarkerEdgeColor','w','LineWidth',5,'SizeData',size_data);
  s2 = scatter3(a,P(:,1),P(:,2),2*offset + P(:,3),'o', ...
    'MarkerFaceColor',c,'MarkerEdgeColor','k','LineWidth',2,'SizeData',size_data);
  s = [s1;s2];
end
