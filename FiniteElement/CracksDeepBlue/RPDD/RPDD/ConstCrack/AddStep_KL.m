function Add=AddStep_KL(Xmarg, X_R)

Add=[];
if iscell(Xmarg)
    for i=1:length(Xmarg)
        xmarg=Xmarg{i};
        SZ=length(xmarg(1,:));
        for ii=1:SZ
            Epsilon_X=KLError( X_R,xmarg(:,ii));
            MaxProbBin=max(Epsilon_X);
            if MaxProbBin==0
                continue
            end
            IAdd=find(Epsilon_X==MaxProbBin);  
            Add=[Add IAdd'];
        end
    end
else
    for i=1:length(Xmarg(1,:))
        Epsilon_X=KLError(X_R,Xmarg(:,i));
        MaxProbBin=max(Epsilon_X);
        if MaxProbBin==0
            continue
        end
        IAdd=find(Epsilon_X==MaxProbBin);
        Add=[Add IAdd'];
    end
    
end

