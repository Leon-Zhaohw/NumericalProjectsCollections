function [Omega,N0,N1,N2,outside_region_of_interest ] = ...
  layers_from_handle( ...
    vertex_count, ...
    F, ...
    H)
  % Given a Face index list and list of vertex indices in a handle (can be
  % disjoint), and number of vertices
  % this function finds the ring or vertices 0, 1, and 2
  % edges into the handle 
  %
  % Input:
  %  vertex_count: number of vertices
  %  F: face index list
  %  H: vertex indices list of vertices in handle
  %
  % Ouput:
  %  Omega: vertex indices list of vertices outside of handle (interior)
  %  N0: vertex indices list of vertices on ring 0 edges into handle
  %  N1: vertex indices list of vertices on ring 1 edges into handle
  %  N2: vertex indices list of vertices on ring 2 edges into handle
  %  outside_region_of_interest: vertex indices list of vertices too far inside
  %                              handle
  %  Omega_1: vertex indices list of vertices one edge into the interior (Omega)
  %
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, NYU
  %

  % omega is simply the whole list minus the handle
  Omega = 1:vertex_count;
  Omega = Omega(~ismember(Omega,H));


  % must have at least vertex in handle and one in Omega
  % interior_faces = F(...
  %   (ismember(F(:,1),H)     & ismember(F(:,2),Omega) & ismember(F(:,3),Omega)) | ...
  %   (ismember(F(:,1),Omega) & ismember(F(:,2),H)     & ismember(F(:,3),Omega)) | ...
  %   (ismember(F(:,1),Omega) & ismember(F(:,2),Omega) & ismember(F(:,3),H)) | ...
  %   (ismember(F(:,1),H)     & ismember(F(:,2),H)     & ismember(F(:,3),Omega)) | ...
  %   (ismember(F(:,1),Omega) & ismember(F(:,2),H)     & ismember(F(:,3),H)) | ...
  %   (ismember(F(:,1),H)     & ismember(F(:,2),Omega) & ismember(F(:,3),H)) ...
  %   ,:);
  % Omega_1 = intersect(Omega,interior_faces(:));

  % Two sets that will grow and shrink as we discover the next layer
  shrinking_handle = H;
  growing_interior = Omega;

  interior_faces = intersect(F,F.*ismember(F,growing_interior),'rows');
  handle_faces = intersect(F,F.*~ismember(F,growing_interior),'rows');
  H0 = setdiff(F,union(handle_faces,interior_faces,'rows'),'rows');
  N0 = intersect(shrinking_handle,reshape(H0,1,size(H0,1)*size(H0,2)));
  growing_interior = [growing_interior N0];
  shrinking_handle = setdiff(shrinking_handle,N0);

  interior_faces = intersect(F,F.*ismember(F,growing_interior),'rows');
  handle_faces = intersect(F,F.*~ismember(F,growing_interior),'rows');
  H1 = setdiff(F,union(handle_faces,interior_faces,'rows'),'rows');
  N1 = intersect(shrinking_handle,reshape(H1,1,size(H1,1)*size(H1,2)));
  growing_interior = [growing_interior N1];
  shrinking_handle = setdiff(shrinking_handle,N1);

  interior_faces = intersect(F,F.*ismember(F,growing_interior),'rows');
  handle_faces = intersect(F,F.*~ismember(F,growing_interior),'rows');
  H2 = setdiff(F,union(handle_faces,interior_faces,'rows'),'rows');
  N2 = intersect(shrinking_handle,reshape(H2,1,size(H2,1)*size(H2,2)));
  growing_interior = [growing_interior N2];
  shrinking_handle = setdiff(shrinking_handle,N2);

  region_of_interest = growing_interior;
  outside_region_of_interest = shrinking_handle;

end
