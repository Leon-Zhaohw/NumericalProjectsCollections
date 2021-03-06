function [YAdd_Predict YE_Predict atAdd_Predict atE_Predict]=AddStep_Predict(atmarg_Jtree,Ymarg_Jtree,at_CPD,T,Inspect)
YAdd_Predict=[];
atAdd_Predict=[];
Dat=cell(1,T-1);
Lat=length(atmarg_Jtree(:,1));
LY=length(Ymarg_Jtree(:,1));
YE_Predict=zeros(LY,1);
atE_Predict=zeros(Lat,1);
for t=max(Inspect(:,1)):T-1
    Dat{1,t}=zeros(Lat-1,LY);
    for i=1:Lat-1
        for j=1:LY
            Dat{1,t}(i,j)=Dat{1,t}(i,j)+sum(atmarg_Jtree(i,t)*Ymarg_Jtree(j,t+1)*at_CPD(i,j,i+1:end));
            atE_Predict(i)=atE_Predict(i)+sum(atmarg_Jtree(i,t)*Ymarg_Jtree(j,t+1)*at_CPD(i,j,i+1:end));
            YE_Predict(j)=YE_Predict(j)+sum(atmarg_Jtree(i,t)*Ymarg_Jtree(j,t+1)*at_CPD(i,j,i+1:end));
        end
    end
    %error is from the maximum probability of at interval that move
    %forward
    MaxProb(t)=max(max(Dat{1,t}));  
    [Imax,Jmax]=find(Dat{1,t}==MaxProb(t)); 
    atAdd_Predict=[atAdd_Predict Imax'];
    YAdd_Predict=[YAdd_Predict Jmax'];
end
atE_Predict=atE_Predict/sum(atE_Predict);
YE_Predict=YE_Predict/sum(atE_Predict);