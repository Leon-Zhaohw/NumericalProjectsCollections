function MergInterval=MergStep(marg_Jtree)

MergInterval=[];
for t=1:length(marg_Jtree(1,:))
    IMerg=find(marg_Jtree(:,t)==0);
    if t==1
        MergInterval=IMerg';
    else
        MergInterval=intersect(MergInterval,IMerg');
    end
end



%function MergInterval=MergStep(marg2_Jtree,Esmooth,Epredict,MergNum)
% for t=1:length(marg2_Jtree)
%     marg2_t=marg2_Jtree{t};
%     if isempty(marg2_t)==0
%         for tt=1:length(marg2_t(1,:))
%              MinProbBin(tt)=min(marg2_t(:,tt));
%              IMerg=find(marg2_t(:,tt)==MinProbBin(tt));
%              if tt==1 && t==1
%                  MergInterval=IMerg';
%              else
%                  MergInterval=intersect(MergInterval,IMerg');
%              end
%         end
%     end
% end
% 
% Error=Esmooth+Epredict;
% Error(end)=inf;
% ErrorComb=Error(1)+Error(2);
% Error(1)=ErrorComb*(1+1E-9);
% Error(2)=ErrorComb*(1-1E-9);
% MergError=Error'
% MergInterval2=[];
% Error(3)=inf;
% Iter=0;    
% 
% 
% while min(Error(6:end))<max(Error(6:end-1))/50 && Iter<MergNum
%     Iter=Iter+1;
%     Merg=find(Error==min(Error(6:end)));
%     if Merg==1
%         Merg=Merg+1;
%     end
%     if isempty(intersect(MergInterval2,Merg+1))&& Error(Merg-1)>Error(Merg+1)
%         Merg=Merg+1;
%     elseif isempty(intersect(MergInterval2,Merg))&& Error(Merg-1)<Error(Merg+1)
%         Merg=Merg;
%     else
%         break
%     end
%     ErrorComb=Error(Merg)+Error(Merg-1);
%     Error(Merg)=ErrorComb*(1+1E-9);
%     Error(Merg-1)=ErrorComb*(1-1E-9);
%     MergInterval2=[MergInterval2 Merg];
% end
% 
% Iter
% MergInterval=[MergInterval MergInterval2];