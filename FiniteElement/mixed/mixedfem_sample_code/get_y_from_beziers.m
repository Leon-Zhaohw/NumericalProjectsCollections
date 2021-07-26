function [y] = get_y_from_beziers(F,V,L,K,BZ1,BZ2)
  % y = Δx =  ∂²x/∂t² + ∂²x/∂n²
  % compute ∂²x/∂n² from quadratic bezier curves (BZ1, BZ2) specified at points
  % in V indexed in L.
  % compute ∂²x/∂t² by computing derivitives of bezier curve formed by a point
  % and its two neighbors: L must be a proper ring.
  %
  % Input:
  %   F: #x3 face list
  %   V: Nx3 vertex list
  %   L: Lx3 index list of vertices in layer (RING)
  %   K: Kx3 index list of vertices in next (outer) layer (RING)
  %   BZ1: Nx3 first vector of bezier curves only L vertices used
  %   BZ2: Nx3 second vector of bezier curves only L vertices used
  %
  % Output:
  %   y: Lx3 values of Δx for all vertices in L 
  %
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, NYU
  %

  number_of_vertices = size(V,1);
  LE = two_rings_to_edge_loops(F,L,K);


  %
  %         p1------p2
  %        /  
  %       /
  %     p0
  %
  % 2nd derivative of quadratic bezier curve at p1 is:
  %   2(p0-p1) + 2(p2-p1)
  %

  d2xdt2 = zeros(number_of_vertices,3);
  d2xdn2 = zeros(number_of_vertices,3);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SCALING TERM: 1/∂n² and 1/∂t²
  dt2 = zeros(number_of_vertices,1);
  dt2(LE(:,1)) = sum((V(LE(:,1),:)-V(LE(:,2),:)).^2,2);
  dt2(LE(:,2)) = dt2(LE(:,2)) + sum((V(LE(:,1),:)-V(LE(:,2),:)).^2,2);
  [valid_L,m,n] = unique(LE(:));
  counts = accumarray(n(:),1);
  % average
  dt2(valid_L) = dt2(valid_L)./counts;

  % also disclud (zero-out) places wehere ∂²x/∂n² is not defined
  valid_L = valid_L(sum(BZ1(valid_L,:).^2,2)~=0);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % tip to tail, tail to tip
  i=[LE(:,1);LE(:,2)]; 
  j=[LE(:,2);LE(:,1)];

  for coord_index = 1:3
    value = ...
      2*[V(LE(:,1),coord_index)-V(LE(:,2),coord_index);...
      V(LE(:,2),coord_index)-V(LE(:,1),coord_index)];
    % sparse takes care of adding up duplicate entries
    %S = sparse(i,j,value);
    S = sparse(i,j,value,number_of_vertices,number_of_vertices);
    % sum up columns
    d2xdt2(:,coord_index) = full(sum(S));


    %
    %
    %      BZ1
    %   p0----->p1
    %            |
    %            | BZ2
    %            ↓
    %            p2
    %                    
    % 2nd derivative of quadratic bezier curve at p1 is:
    %   2(p0-p1) + 2(p2-p1)
    %   2(-BZ1)  + 2(BZ2)
    %d2xdn2(:,coord_index) = -2 * BZ1(:,coord_index) + 2 * BZ2(:,coord_index);
    % don't know why this is working better than above
    d2xdn2(:,coord_index) = 2 * BZ1(:,coord_index) + 2 * BZ2(:,coord_index);

  end

  % Δx =  ∂²x/∂t² + ∂²x/∂n²
  %y = d2xdt2(L,:) + d2xdn2(L,:);
  y = zeros(number_of_vertices,3);
  y(valid_L,:) = d2xdt2(valid_L,:)./repmat(dt2(valid_L),1,3) + d2xdn2(valid_L,:)./repmat(dt2(valid_L),1,3);
  %error(num2str(max(dt2(valid_L,:).^-1)))
  y = y(L,:);
  %error(num2str(max(1.0/dt2)))
  % somehow, they're backwards.....
  %y = -d2xdt2(L,:) - d2xdn2(L,:);
end


