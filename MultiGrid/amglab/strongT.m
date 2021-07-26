function Si = strongT(i,V,S,level);
% Determine the transposed set of strong connections [Ruge/Stüben].
% SCPT(I,V,S) returns the transposed (and sorted) set of strong connections
% of node I from the remaining set of free nodes V 
% with respect to the already calculated 
% array of strong connection sets S.S

sk_set=[]; % initialize the transposed strong set to be empty
for k=V % consider each point in the domain
    if k~=i % do not consider connections to one's self...
       bb = ismember(i,S(k).S); % find out if point I is a member of the strongly connected set of k
       if bb(1) % if it is...
           sk_set=union(sk_set,k); % add K to the strongly connected set of I
       end
    end
end
Si = sort(sk_set); % sort and return the transposed strongly connected set of I