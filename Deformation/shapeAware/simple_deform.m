function varargout = simple_deform(varargin)
  % SIMPLE_DEFORM Deform a mesh (V,F) using linear blend skinning or dual quaternion
  % skinning, by applying transformations at a set of control points C and
  % propogating the deformation to the mesh via correspondence weights W.
  %
  % simple_deform(V,F,C,W,...)
  %
  % Inputs:
  %  V  #V by 2 list of vertex positions
  %  F  #F by 3 list of face indices
  %  C  #C by 2 list of control point positions
  %  W  weights, #vertices by #handles matrix of weights
  %  OPTIONS:
  %    'PointHandles' followed by a list of indices into C, for the point
  %      handles. Default value is 1:size(C,1)
  %    'BoneEdges' followed by a #BE by 2 list of indices into C for the bone
  %      handles. Default value is []
  %    'CageEdges' followed by a #CE by 2 list of indices into P specifying 
  %      cage edges for point handles (for display purposes only)
  %    'ShowWeightVisualization' puts a subplot to the right with the selected
  %      handle's weights in a plot using color to indicate value, default is 
  %      off
  %    'InterpMode' could be 'LBS' for Linear Blend Skinning or 'DQLBS' for
  %      Dual Quaternion Linear Blend Skinning, default is 'LBS'
  %    'AutoDOF' followed by a method to compute automatic degrees of freedom
  %      for point handles, possible options are:
  %       'none'
  %       'pseudoedges'
  %    'PseudoEdges' followed by a #PE by 2 list of pseudo edges for computing
  %      automatic rotations at point handles
  % Output:
  %   gid  index into global variable g_Deform, which gives access to plot
  %     handles and input variables
  %
  % Global Output:
  %   struct array g_Deform:
  %     g_Deform(gid).R  pose rotations of point handles
  %     g_Deform(gid).new_C  pose positions of control vertices
  %     g_Deform(gid).update_positions  function handle to update positions
  %       based on fields R and new_C
  %     g_Deform(gid).tsh  plot handle to main triangle mesh plot
  %     g_Deform(gid).wvsh  plot handle to weight visualization plot
  % 
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: deform
  %

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % parse input
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % vertex positions of mesh
  V = varargin{1};
  % face indices of mesh
  F = varargin{2};
  % control vertices of skeleton
  C = varargin{3};
  % deformation weights
  W = varargin{4};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set default parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % set default point handles
  P = 1:size(C,1);
  % set default bone edge handles
  BE = [];
  % Be sure that control vertices are in 2D
  if(size(C,2) == 3)
    C = C(:,1:2);
  end
  % Weights used for contours
  CW = W;
  % Weights used for weight visualization
  WVW = W;
  % show weight visualization off to right
  show_weight_visualization = false;
  % use linear blend skinning to deform (other option is DQLBS)
  interp_mode = 'LBS';
  % use automatically computed dof at point handles
  auto_dof = 'none';
  % pseudo-edges
  PE = [];
  % cage edges (for display only)
  CE = [];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Parse additional options
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii = 5;
  while(ii <= size(varargin,2))
    if(strcmp(varargin{ii},'ShowWeightVisualization'))
      show_weight_visualization = true;
    elseif(strcmp(varargin{ii},'InterpMode'))
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      if( strcmp(varargin{ii},'LBS')) 
        interp_mode = 'LBS';
      elseif( strcmp(varargin{ii},'DQLBS')) 
        interp_mode = 'DQLBS';
      else
        error('InterpMode must be either LBS or DQLBS');
      end
    elseif(strcmp(varargin{ii},'BoneEdges'))
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      BE = varargin{ii};
    elseif(strcmp(varargin{ii},'CageEdges'))
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      CE = varargin{ii};
    elseif(strcmp(varargin{ii},'PointHandles'))
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      P = varargin{ii};
    elseif(strcmp(varargin{ii},'AutoDOF'))
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      auto_dof = varargin{ii};
    elseif(strcmp(varargin{ii},'PseudoEdges'))
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      PE = varargin{ii};
    else
      error('Bad parameter');
    end
    ii = ii + 1;
  end

  % number of point handles
  np = numel(P);
  % number of bone handles
  nb = size(BE,1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Prepare output
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  global g_Deform;
  gid = numel(g_Deform)+1;
  % keep track of control positions at mouse down
  g_Deform(gid).new_C = [];
  % keep track of rotations stored at each control point, for 2D this is a m
  % by 1 list of angles
  g_Deform(gid).R = zeros(np,1);
  g_Deform(gid).update_positions = @update_positions;
  % set output if it exists
  if(nargout == 0)
  elseif(nargout == 1)
    varargout{1} = gid;
  else
    error('Wrong number of output arguements');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set up plots
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Set up figure with enough subplots for any additional visualizations
  % clear current figure
  clf
  number_of_subplots = 1;
  current_subplot = 1;
  if(show_weight_visualization)
    number_of_subplots = number_of_subplots +1;
  end
  if(number_of_subplots>1)
    subplot(1,number_of_subplots,1);
  end

  % plot the original mesh
  g_Deform(gid).tsh = trisurf(F,V(:,1),V(:,2),zeros(size(V,1),1), ...
    'FaceColor','interp');
  % 2D view
  view(2);
  axis equal
  axis manual
  % plot bones
  hold on;
  if(nb > 0)
    % plot thick lines for bones (outline of lines)
    B_plot_outer = plot( ...
      [C(BE(:,1),1) C(BE(:,2),1)]', ...
      [C(BE(:,1),2) C(BE(:,2),2)]', ...
      '-k', ...
      'LineWidth',5);
    % plot thin lines for bones (innerline of lines)
    B_plot_inner = plot( ...
      [C(BE(:,1),1) C(BE(:,2),1)]', ...
      [C(BE(:,1),2) C(BE(:,2),2)]', ...
      '-b', ...
      'LineWidth',2);
  end
  if(size(PE,1) > 0)
    % plot lines for pseudo edges
    PE_plot = plot( ...
      [C(P(PE(:,1)),1) C(P(PE(:,2)),1)]', ...
      [C(P(PE(:,1)),2) C(P(PE(:,2)),2)]', ...
      '--r', ...
      'LineWidth',5);
  end
  if(size(CE,1) > 0)
    % plot lines for cage edges
    CE_plot_outer = plot( ...
      [C(P(CE(:,1)),1) C(P(CE(:,2)),1)]', ...
      [C(P(CE(:,1)),2) C(P(CE(:,2)),2)]', ...
      '-k', ...
      'LineWidth',5);
    CE_plot_inner = plot( ...
      [C(P(CE(:,1)),1) C(P(CE(:,2)),1)]', ...
      [C(P(CE(:,1)),2) C(P(CE(:,2)),2)]', ...
      '-', ...
      'Color', [1 1 0.2], ...
      'LineWidth',2);
  end
  % plot the control points (use 3D plot and fake a depth offset by pushing
  % control points up in z-direction)
  C_plot = scatter3( ...
    C(:,1),C(:,2),0.1+0*C(:,1), ... 
    'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
    'LineWidth',2,'SizeData',100, ...
    'ButtonDownFcn',@oncontrolsdown);
  hold off;

  % set up weight visualization plot
  if(show_weight_visualization)
    current_subplot = current_subplot + 1;
    assert(all(size(CW) == size(W)));
    % subplot for weight visualization
    subplot(1,number_of_subplots,current_subplot);
    hold on;
    % plot the original mesh
    g_Deform(gid).wvsh = ...
      trisurf(F,V(:,1),V(:,2),zeros(size(V,1),1), ...
        'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    view(2);
    axis equal
    axis manual
    hold off;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set up interaction variables
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % keep track of window xmin, xmax, ymin, ymax
  win_min = min([C(:,1:2); V(:,1:2)]);
  win_max = max([C(:,1:2); V(:,1:2)]);
  % keep track of down position
  down_pos = [];
  % keep track of last two drag positions
  drag_pos = [];
  last_drag_pos = [];
  % keep track of mesh vertices at mouse down
  down_V = [];
  % keep track of index of selected control point
  ci = [];
  % type of click ('left','right')
  down_type  = '';

  if(show_weight_visualization)
    fprintf(['\nCLICK a control point to visualize its corresponding ' ...
      'weights on the mesh.\n']);
  end
  fprintf( ...
    ['DRAG a control point to deform the mesh.\n' ...
    'RIGHT CLICK DRAG a control point to rotate point handles.\n\n']);

  return

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Callback functions for keyboard and mouse
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Callback for mouse down on control points
  function oncontrolsdown(src,ev)
    % get current mouse position, and remember old one
    down_pos=get(gca,'currentpoint');
    down_pos=[down_pos(1,1,1) down_pos(1,2,1)];
    last_drag_pos=down_pos;
    drag_pos=down_pos;
    % keep track of control point positions at mouse down
    g_Deform(gid).new_C = [get(C_plot,'XData')' get(C_plot,'YData')'];
    % get index of closest control point
    [minD,ci] =  ...
      min(sum((g_Deform(gid).new_C(:,1:2) - ...
      repmat(down_pos,size(g_Deform(gid).new_C,1),1)).^2,2));
    % keep track of mesh vertices at mouse down
    down_V = get(g_Deform(gid).tsh,'Vertices');
    down_V = down_V(:,1:2);

    % tell window that drag and up events should be handled by controls
    set(gcf,'windowbuttonmotionfcn',@oncontrolsdrag)
    set(gcf,'windowbuttonupfcn',@oncontrolsup)
    set(gcf,'KeyPressFcn',@onkeypress)
    if(strcmp('normal',get(gcf,'SelectionType')))
      % left-click
      down_type = 'left';
    else
      % other (right) click
      down_type = 'right';
    end

    % try to find ci in list of point handles
    [found, iP] = ismember(ci,P);
    if(found)
      % set color of mesh plot to weights of selected
      %set(tsh,'CData',W(:,iP));
      % change weights in weight visualization
      if(show_weight_visualization)
        set(g_Deform(gid).wvsh,'CData',WVW(:,iP));
      end
    end

  end

  % Callback for mouse drag on control points
  function oncontrolsdrag(src,ev)
    % keep last drag position
    last_drag_pos=drag_pos;
    % get current mouse position
    drag_pos=get(gca,'currentpoint');
    drag_pos=[drag_pos(1,1,1) drag_pos(1,2,1)];
    if(strcmp('left',down_type))
      % move selected control point by drag offset
      g_Deform(gid).new_C(ci,:) = ...
        g_Deform(gid).new_C(ci,:) + drag_pos-last_drag_pos;
    else
      [found, iP] = ismember(ci,P);
      if(found)
        g_Deform(gid).R(iP) = ...
          g_Deform(gid).R(iP) + 2*pi*(drag_pos(1)-last_drag_pos(1))/100;
      end
    end
    update_positions();
  end

  function update_positions()
    % update display positions
    set(C_plot,'XData',g_Deform(gid).new_C(:,1));
    set(C_plot,'YData',g_Deform(gid).new_C(:,2));
    if(nb > 0)
      set(B_plot_outer,{'XData'}, num2cell([ ...
        g_Deform(gid).new_C(BE(:,1),1) ...
        g_Deform(gid).new_C(BE(:,2),1)],2));
      set(B_plot_outer,{'YData'}, num2cell([ ...
        g_Deform(gid).new_C(BE(:,1),2) ...
        g_Deform(gid).new_C(BE(:,2),2)],2));
      set(B_plot_inner,{'XData'}, num2cell([ ...
        g_Deform(gid).new_C(BE(:,1),1) ...
        g_Deform(gid).new_C(BE(:,2),1)],2));
      set(B_plot_inner,{'YData'}, num2cell([ ...
        g_Deform(gid).new_C(BE(:,1),2) ...
        g_Deform(gid).new_C(BE(:,2),2)],2));
    end
    % update pseudo edge plots
    if(size(PE,1)>0)
      set(PE_plot,{'XData'}, num2cell([ ...
        g_Deform(gid).new_C(P(PE(:,1)),1) ...
        g_Deform(gid).new_C(P(PE(:,2)),1)],2));
      set(PE_plot,{'YData'}, num2cell([ ...
        g_Deform(gid).new_C(P(PE(:,1)),2) ...
        g_Deform(gid).new_C(P(PE(:,2)),2)],2));
    end
    % update cage edge plots
    if(size(CE,1)>0)
      set(CE_plot_outer,{'XData'}, num2cell([ ...
        g_Deform(gid).new_C(P(CE(:,1)),1) ...
        g_Deform(gid).new_C(P(CE(:,2)),1)],2));
      set(CE_plot_outer,{'YData'}, num2cell([ ...
        g_Deform(gid).new_C(P(CE(:,1)),2) ...
        g_Deform(gid).new_C(P(CE(:,2)),2)],2));
      set(CE_plot_inner,{'XData'}, num2cell([ ...
        g_Deform(gid).new_C(P(CE(:,1)),1) ...
        g_Deform(gid).new_C(P(CE(:,2)),1)],2));
      set(CE_plot_inner,{'YData'}, num2cell([ ...
        g_Deform(gid).new_C(P(CE(:,1)),2) ...
        g_Deform(gid).new_C(P(CE(:,2)),2)],2));
    end
    % update mesh positions

    % USING LINEAR BLEND SKINNING
    if(strcmp(interp_mode,'LBS'))
      % get transformations stored at each point and bone handle
      TR = ...
        skinning_transformations(C,P,BE,g_Deform(gid).new_C,g_Deform(gid).R);
      if strcmp(auto_dof,'pseudoedges')
        TR(:,:,1:np) = pseudoedge_dof( ...
          C(P,:),PE,g_Deform(gid).new_C(P,:)-C(P,:), ...
          axisangle2quat(repmat([0 0 1],np,1),g_Deform(gid).R));
      end

      % linear blend skinning
      [new_V] = lbs(V(:,1:2),TR,W);
    elseif(strcmp(interp_mode,'DQLBS'))

      % USING DUAL QUATERNION SKINNING
      % convert angles around z-axis to quaternions
      [T,AX,AN,Sm,O] = ...
        skinning_transformations(C,P,BE,g_Deform(gid).new_C,g_Deform(gid).R);
      % number of handles
      m = numel(P)+size(BE,1);
      % dimension (2 or 3)
      dim = size(C,2);
      % Extract scale
      TR = zeros(dim,dim+1,m);
      TR(1:dim,1:dim,:) = Sm;
      Sm = reshape(Sm,[dim dim*m])';
      TR(1:dim,dim+1,:) = permute(O-stacktimes(Sm,O),[2 3 1]);
      % Perform scale as linear blend skinning, before translations and
      % rotations
      [new_V] = lbs(V(:,1:2),TR,W);

      Q = axisangle2quat(AX,AN);
      % quattrans2udq expect 3D translations, so pad with zeros
      T = [T zeros(size(T,1),1)];
      % convert quaternions and translations into dualquaternions
      DQ = quattrans2udq(Q,T);

      new_V = dualquatlbs(new_V,DQ,W);
    end

    % update mesh positions
    set(g_Deform(gid).tsh,'Vertices',new_V(:,1:2));
  end

  % Callback for mouse release of control points
  function oncontrolsup(src,ev)
    % Tell window to handle drag and up events itself
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
    cur_V = get(g_Deform(gid).tsh,'Vertices');
    cur_V = cur_V(:,1:2);

    % scale window to fit
    win_min = min([win_min; cur_V]);
    win_max = max([win_max; cur_V]);
    axis(reshape([win_min;win_max],1,2*size(cur_V,2)))
  end

  function onkeypress(src,ev)
    if(strcmp(ev.Character,'r'))
      g_Deform(gid).new_C = C;
      g_Deform(gid).R = zeros(np,1);
      update_positions();
    elseif(strcmp(ev.Character,'u'))
      update_positions();
    end
  end

end
