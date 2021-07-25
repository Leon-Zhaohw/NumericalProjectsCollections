function [Add_Smooth Error_Smooth]=AddStep_Smooth(marg_Jtree)
Add_Smooth=[];
L=length(marg_Jtree{1}(:,1));
Error_Smooth=zeros(L,1);
for i=1:length(marg_Jtree)
    marg_t=marg_Jtree{i};
    T=length(marg_Jtree{1}(1,:));
    for ii=1:T
        Error_Smooth=Error_Smooth+marg_t(:,ii);
         MaxProbBin=max(marg_t(:,ii));
         IAdd=find(marg_t(:,ii)==MaxProbBin);
         Add_Smooth=[Add_Smooth IAdd'];
    end
end
Error_Smooth=Error_Smooth/sum(Error_Smooth);

