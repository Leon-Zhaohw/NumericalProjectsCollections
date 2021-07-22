function Volkraft = f(u,zeitpkt);
% Volumenkraft des Problems
% u Matrix, Punkte zeilenweise

 if zeitpkt<0
   Volkraft = zeros(size(u,1),2);
 else
   Volkraft = (1-zeitpkt)*zeitpkt*zeros(size(u,1),2);
 end
