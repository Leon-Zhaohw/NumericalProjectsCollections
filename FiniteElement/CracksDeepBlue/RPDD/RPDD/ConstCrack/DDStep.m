function Range=DDStep(AddInterval, MergInterval,Range)        
%Dynamic Discretization operation which includes add and merge
% Range is the discretization for both input and output
AddInterval=sort(AddInterval,'descend');
MergInterval=sort(MergInterval,'descend');
L=length(Range);
AddFlag=-1;
for i=AddInterval
    if i<=L && i>=1
        AddFlag=1;
    end
end
MergFlag=-1;
if isempty(MergInterval)==0
    for i=MergInterval
        if i<=L && i>=1
            MergFlag=1;
        end
    end
end
if isempty(intersect(AddInterval,MergInterval))==0
    % if Add Intervals and Merge Intervals have intersection then do
    % nothing about that interval
    SharedIntervals=intersect(AddInterval,MergInterval);
    for sharedinterval=SharedIntervals
        Iadd=find(AddInterval==sharedinterval);
        AddInterval(Iadd)=[];
        Imerg=find(MergInterval==sharedinterval);
        MergInterval(Imerg)=[];
    end
end

if AddFlag==1 && MergFlag==1
    % If need both add and merge operation
    Index=unique([MergInterval AddInterval]);
    Index=sort(Index,'descend');
    for I=Index
        if isempty(find(AddInterval==I))==0
            %Add Operation
            Range=[Range(1:I) (Range(I)+Range(I+1))/2 Range(I+1:end)];
        else
            %Merg Operation
            Range(I)=[]; 
        end
    end
elseif AddFlag==1 && MergFlag~=1
    %iIf just add operation
    for I=AddInterval
        Range=[Range(1:I) (Range(I)+Range(I+1))/2 Range(I+1:end)];           
    end
elseif AddFlag~=1 && MergFlag==1
    % If just merge operation
    for I=MergInterval
        Range(I)=[];
    end
end

    
    
    
    
    
    
    
    
    