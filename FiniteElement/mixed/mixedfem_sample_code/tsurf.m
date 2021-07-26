function tsurf(F,V,vertex_indices,face_indices)
  % trisurf wrapper
  % F: list of faces #F x 3
  % V: vertex positiosn #V x 3 or #V x 2
  % vertex_indices: show vertex indices on plot
  %                 0 -> off
  %                 1 -> text and grey background
  %                >1 -> text
  % face_indices: show face indices on plot
  %                 0 -> off
  %                 1 -> text and grey background
  %                >1 -> text
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, NYU
  %

  if(~exist('vertex_indices'))
    % off by default
    vertex_indices = 0;
  end
  if(~exist('face_indices'))
    % off by default
    face_indices = 0;
  end

  is_2d = false;
  if(size(V,2)==2 || (size(V,2) ==3 && sum(abs(V(:,3))) == 0))
    V = [V(:,1) V(:,2) 0*V(:,1)];
    is_2d = true;
  elseif(size(V,2)>3 || size(V,2)<2 ) 
    error('V must be #V x 3 or #V x 2');
    return;
  end

  trisurf(F,V(:,1),V(:,2),V(:,3));

  % if 2d then set to view (x,y) plane
  if( is_2d )
    view(2);
  end

  if(face_indices==1)
    FC = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))./3;
    text(FC(:,1),FC(:,2),FC(:,3),num2str((1:size(F,1))'),'BackgroundColor',[.8 .8 .8]);
  elseif(face_indices)
    FC = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))./3;
    text(FC(:,1),FC(:,2),FC(:,3),num2str((1:size(F,1))'));
  end

  if(vertex_indices==1)
    text(V(:,1),V(:,2),V(:,3),num2str((1:size(V,1))'),'BackgroundColor',[.8 .8 .8]);
  elseif(vertex_indices)
    text(V(:,1),V(:,2),V(:,3),num2str((1:size(V,1))'));
  end
  % uncomment these to switch to a better 3d surface viewing mode
  %axis equal; axis vis3d;

end
