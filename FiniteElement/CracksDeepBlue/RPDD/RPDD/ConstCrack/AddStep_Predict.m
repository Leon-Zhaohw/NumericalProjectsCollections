function  [atAdd_Predict daAdd_Predict atE_Predict daE_Predict]=AddStep_Predict(damarg_Jtree,atmarg_Jtree,at_CPD,T,Evidence)
atAdd_Predict=[];
daAdd_Predict=[];
Dat=cell(1,T-1);
Lda=length(damarg_Jtree(:,1));
Lat=length(atmarg_Jtree(:,1));
atE_Predict=zeros(Lat,1);
daE_Predict=zeros(Lda,1);
%for t=max(Evidence(:,1))+1:T-1
for t=max(Evidence(:,1)):T-1
    Dat{1,t}=zeros(Lat-1,Lda);
    for i=1:Lat-1
        for j=1:Lda
            Dat{1,t}(i,j)=Dat{1,t}(i,j)+sum(atmarg_Jtree(i,t)*damarg_Jtree(j,t+1)*at_CPD(i,j,i+1:end));
            atE_Predict(i)=atE_Predict(i)+sum(atmarg_Jtree(i,t)*damarg_Jtree(j,t+1)*at_CPD(i,j,i+1:end));
            daE_Predict(j)=daE_Predict(j)+sum(atmarg_Jtree(i,t)*damarg_Jtree(j,t+1)*at_CPD(i,j,i+1:end));
        end
    end
    %error is from the maximum probability of at interval that move
    %forward
    MaxProb(t)=max(max(Dat{1,t}));  
    [Imax,Jmax]=find(Dat{1,t}==MaxProb(t)); 
    atAdd_Predict=[atAdd_Predict Imax'];
    daAdd_Predict=[daAdd_Predict Jmax'];
    
end
atE_Predict=atE_Predict/sum(atE_Predict);
daE_Predict=daE_Predict/sum(daE_Predict);