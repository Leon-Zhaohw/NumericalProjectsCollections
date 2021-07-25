function Range=DDalgorithm(marg_ev,Range,IndexProtect)
% Dynamic Discretisation algorithm by Neil& Fenton
Range_DD=Range;
%assume a interval size for last interval
Range_DD(end)=Range(end-1)+10;
L=length(marg_ev);
marg_ev
if length(marg_ev)>1
    for i=1:L
        if i==1
            fmin=0;
            fmax=marg_ev(2);
            fbar=((Range_DD(2)-Range_DD(1))*marg_ev(1)+(Range_DD(3)-Range_DD(2))*fmax)/2;
        elseif i==L
            fmax=marg_ev(i-1);
            fmin=0;
            fbar=((Range_DD(L)-Range_DD(L-1))*fmax+(Range_DD(L+1)-Range_DD(L))*marg_ev(L))/2;
        else
            f1=marg_ev(i-1);
            f2=marg_ev(i);
            f3=marg_ev(i+1);
            fbar=( (Range_DD(i)-Range_DD(i-1))*f1 + (Range_DD(i+1)-Range_DD(i))*f2+(Range_DD(i+2)-Range_DD(i+1))*f3)/3;
            fmin=min([f1 f2 f3]);
            fmax=max([f1 f2 f3]);
        end
        Ej(i)=(fmax-fbar)/(fmax-fmin)*fmin*log(fmin/fbar)+(fbar-fmin)/(fmax-fmin)*fmax*log(fmax/fbar)*abs(Range_DD(i+1)-Range_DD(i));
    end
    Ej
    AddInterval=find(Ej==max(Ej));
    MergInterval=[find(isnan(Ej)==1)];% find(Ej==min(Ej))];
    % Merge with the smallest neighbour
    MergDirection=MergInterval;
    for i=1:length(MergInterval)
        if MergInterval(i)==1
            MergDirection(i)=1;
        elseif MergInterval(i)==L
            MergDirection(i)=MergInterval(i);
        elseif Ej(MergInterval(i)+1)<Ej(MergInterval(i)-1)
                MergDirection(i)=MergInterval(i)+1;
        end
    end
    MergDirection=unique(MergDirection);
    I=find(MergDirection==IndexProtect);
    MergDirection(I)=[];
    Range_DD=atDDStep(AddInterval, MergDirection,Range_DD); 
    Range=Range_DD;
    Range(end)=inf;
end




























