function Vals_R=FindInt(RANG,Vals_act)
%Find the interval where the actual value locates.
Vals_R=Vals_act;
for i=1:length(Vals_act)
    if Vals_act{i}~=0
        for j=2:length(RANG)
            if RANG(j)>Vals_act{i}
                Vals_R{i}=j-1;
                break
            end
        end
    end
 end