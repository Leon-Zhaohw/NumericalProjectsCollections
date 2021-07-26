function [toodeep] = RugeStuebenCoarsen(level)
% This is the Ruge Stueben algorithm from the
% McCormick Book "Multigrid Methods."
%
% Derrick Cerwinsky
% University of Wyoming

amg_globals;


if DEBUG == 1
    disp('Starting RugeStuebenCoarsen')
end


if issparse(A(level).matrix) == false    % Test if A(level).matrix is sparse, make sparse if not.
    A(level).matrix = sparse(A(level).matrix);
end


theta = THETA;  % This is set in globals


[m n] = size(A(level).matrix);  % Number of nodes.  The matrix is square, so we do not keep the m.

m = nnz(A(level).matrix);  % Number of non-zero elements in A(level).matrix.

I = zeros( 1, n );  % Initialize the set of nodes.  0 = undeclared, 1 = coarse, 2 = fine.

c = 0;    % Counts coarse nodes.
f = 0;    % Counts fine nodes.

[row, col, a_ij] = find(A(level).matrix);

 cnt = zeros(1,n);  % Count of non-zero elements in each row.
 dsp = zeros(1,n);  % Displacement of elements in cnt.
 
 
 
% This counts the number of non-zero elements in each row.

for i = 1:n
    cnt(i) = nnz(A(level).matrix(i,:));
end

% This finds displacement of elements in cnt.
dsp(1) = 1;
for i = 1:n-1
    dsp(i+1) = dsp(i) + cnt(i);
end


neg_max_val = zeros(n,1);  % This will pick out the max value in each row.


% This will pick out the magnitude of the "largest"
% off diagonal element.  It assumes, however, that
% A(level).matrix is a M matrix.

for i=1:n
   neg_max_val(i) = -min(A(level).matrix(i,:)); 
end


% lambda is set to -1 for unassined nodes.
% It is set to -2 for course nodes.

lambda = ones(1, n)*(-1);



 s_cnt = zeros(1,n);  % Count of non-zero elements in each row.
 s_dsp = zeros(1,n+1);  % Displacement of elements in cnt.
 


S = zeros(n,2);

snnz=0;


% This will fill out the matrix S.

k = 1;

for i = 1:n
   for j = dsp(i):(dsp(i) + cnt(i) - 1) 
       if abs(a_ij(k)) >= theta * neg_max_val(i) && row(k) ~= col(k)
           snnz = snnz + 1;
           S(snnz,:) = [row(k),col(k)];
           if lambda(row(k)) == -1
               lambda(row(k))=0;
           end
           lambda(row(k)) = lambda(row(k)) + 1;
           s_cnt(row(k)) = s_cnt(row(k)) + 1;
       end
       k = k + 1;
   end
end

S = sortrows(S);

% This finds displacement of elements in S.
s_dsp(1) = 1;
for i = 1:n
    if i < n
        s_dsp(i+1) = s_dsp(i) + s_cnt(i);
    else
        s_dsp(i+1) = nnz(S)/2;
    end
end



St = zeros(n,2);

stnnz=0;

st_cnt = zeros(1,n);  % Count of non-zero elements in each row.
st_dsp = zeros(1,n+1);  % Displacement of elements in cnt.



% This will fill out the matrix St.

k = 1;

for i = 1:n
   for j = dsp(i):(dsp(i) + cnt(i) - 1) 
       if abs(a_ij(k)) >= theta * neg_max_val(i) && row(k) ~= col(k)
           stnnz = stnnz + 1;
           St(stnnz,:) = [ col(k),row(k)];
           st_cnt(row(k)) = st_cnt(row(k)) + 1;
       end
       k = k + 1;
   end
end

St = sortrows(St);

 

% This finds displacement of elements in st_cnt.
st_dsp(1) = 1;
for i = 1:n
    if i < n
        st_dsp(i+1) = st_dsp(i) + st_cnt(i);
    else
        st_dsp(i+1) = nnz(St)/2;
    end
end


k = 0;

while k < n
  
  [a i] = max(lambda); % a is the value of the max, i is the index of the max.
  lambda(i) = -2;
  
  I(i) = 1;
  k = k + 1;
  c = c + 1;
  
% Here we start to choose coarse nodes and make all nodes conneded to
% coarse nodes fine nodes.


  for j = st_dsp(i):st_dsp(i+1)-1
      if I(St(j,2)) == 0
          I(St(j,2)) = 2;
          lambda(St(j,2)) = -2;
          f = f + 1;
          k = k + 1;
          for l = s_dsp(St(j,2)):s_dsp(St(j,2)+1) - 1
              if I(S(l,2)) == 0
                  if lambda(S(l,2)) == -1
                      lambda(S(l,2)) = 0;
                  end
              lambda(S(l,2)) = lambda(S(l,2)) + 1;
              end
          end
% This code was moved here to make sure it will run
% However, this code may be redundent.
%             if lambda(S(j,2)) == -1;
%                 lambda(S(j,2)) = 0;
%             end
%             lambda(S(j,2)) = lambda(S(j,2)) + 1;
% End moved code
      end
  end

% This was commented out and moved up.
%    for j = s_dsp(i):s_dsp(i+1)-1
%       if I(S(j,2)) == 0
%           if lambda(S(j,2)) == -1;
%               lambda(S(j,2)) = 0;
%           end
%           lambda(S(j,2)) = lambda(S(j,2)) + 1;
%       end
%    end
end
  
      


% This completes the initial pass.  Next we need to check for F - F
% connections which do not depend on a common C node.  For this we will look at the matrix
% A(level).matrix.    



for i = 1:n
    if I(i) == 2
        flag1 = 1;
        flag2 = 0;
		for j = dsp(i):dsp(i + 1) - 1
            if I(row(j)) == 2 && row(j) ~= col(j)
                flag1=0;
                for k = 1:n
                    if (I(k) == 1 && A(level).matrix(i,k)~=0) && (A(level).matrix(col(j),k)~=0)
                        flag2 = 1;
                    end
                end
			end
		end
		if flag1 == 0 && flag2 == 0
			I(i) = 1;
			f = f - 1;
			c = c + 1;
		end
	end
end

% Now we need to make the interpolation matrix.
% To do this, we need the omega matrix.
% This will have the form omega = - (a_ij + B)/ (a_ii + C)
% Where B = sum ( a_im * a_mj / D)
% C = sum a_in
% D = sum a_mk


Z = zeros(f,c);
x=0;


for i = 1:n
	if I(i) == 2
        x = x + 1;
        y = 0;
        for j = 1:n
            if I(j) == 1
                y = y + 1;
                B = 0;
                % *******************************************
                for m = s_dsp(i):s_dsp(i+1)-1
                    D = 0;
                    % ***************************************
                    
                    
                    %if S(i,m) == 1
                        for k = 1:n
                            if I(k) == 1
                                D = D + A(level).matrix(col(m),k);
                            end
                        end
                    %end
                    
                    B = B + (A(level).matrix(i,row(m)) * A(level).matrix(col(m),j))/ D;
                end
                % *******************************************
                C = 0;
                for n_1 = 1:n
                    if I(n_1) == 2
                        C = C + A(level).matrix(i, n_1);
                    end
                end
		
                Z(x,y) = - (A(level).matrix(i,j) + B) / ( A(level).matrix(i,i) + C );
            end
        end
	end
end

% Now to make the prolongation matrix

%P = zeros(n, c);
P = [0];
P = sparse(P);
P(n,c) = 0;


index_c = 0;
index_f = 0;

for i = 1:n
	if I(i) == 1
		index_c = index_c + 1;
		P(i,index_c) = 1;
	else
		index_f = index_f + 1;
		P(i,:) = Z(index_f , :);
	end
end



% Now to make and return the coarse matrix.

A(level + 1).matrix = sparse(P' * A(level).matrix * P );
W(level).Iweight = P;
W(level).Rweight = P';


if DEBUG == 1
    disp('Exit RugeStuebenCoarsen')
end


toodeep = 0;
