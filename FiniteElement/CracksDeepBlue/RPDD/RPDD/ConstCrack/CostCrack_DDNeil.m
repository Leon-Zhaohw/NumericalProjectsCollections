function    [a0marg_DDNeil, damarg_DDNeil, atmarg_DDNeil, Mmarg_DDNeil, da_R, at_R, M_R]...
        =CostCrack_DDNeil(mu_a0, sigma_a0, mu_da, sigma_da, sigma_error,Evidence, atc, T, da_R,at_R,M_R,a0_R,Lda_Targ,Lat_Targ, LM_Targ)

NumOfSplit=5;
NumOfMerge=0;

a0Add=[];
daAdd=cell(1,T);
atAdd=cell(1,T);
MAdd=cell(1,T);
a0Merg=[];
daMerg=cell(1,T);
atMerg=cell(1,T);
MMerg=cell(1,T);


TotalIter=29;
for Iter=1:TotalIter
    [a0_CPD,da_CPD, at_CPD,M_CPD]=CPDgen_SBN(a0_R,da_R,at_R,M_R,mu_a0,sigma_a0,mu_da, sigma_da,sigma_error,atc);
    ev=cell(1,T*3+1);
    SizeEV=size(Evidence);
    for i=1:SizeEV(1)
        LM=length(M_CPD{Evidence(i,1)}(1,:));
        for j=1:LM
            if M_R{Evidence(i,1)}(j+1)>=Evidence(i,2)
                ev{1,Evidence(i,1)*3+1}=j;
                break
            end
        end
    end  
    [SBNengine,a0marg_DDNeil,damarg_DDNeil, atmarg_DDNeil, Mmarg_DDNeil]=ConstCrack_SBNModel(T,ev,a0_CPD,da_CPD,at_CPD,M_CPD);
    a0_Rmid=[(a0_R(2:end-1)+a0_R(1:end-2))/2 a0_R(end-1)];
    mu_a0_DDNeil=sum(a0_Rmid.*a0marg_DDNeil');
    sigma2_a0_DDNeil=sqrt(sum(a0_Rmid.^2.*a0marg_DDNeil')-mu_a0_DDNeil^2);
    for t=1:T
        da_Rmid=[(da_R{t}(2:end-1)+da_R{t}(1:end-2))/2 da_R{t}(end-1)];
        mu_da_DDNeil(t)=sum(da_Rmid.*damarg_DDNeil{t}');
        sigma2_da_DDNeil(t)= (sum(da_Rmid.^2.*damarg_DDNeil{t}')-mu_da_DDNeil(t)^2);
        at_Rmid=[(at_R{t}(2:end-1)+at_R{t}(1:end-2))/2 at_R{t}(end-1)];
        mu_at_DDNeil(t)=sum(at_Rmid.*atmarg_DDNeil{t}');
        sigma2_at_DDNeil(t)=(sum(at_Rmid.^2.*atmarg_DDNeil{t}')-mu_at_DDNeil(t)^2);
        M_Rmid=[(M_R{t}(2:end-1)+M_R{t}(1:end-2))/2 M_R{t}(end-1)];
        mu_M_DDNeil(t)=sum(M_Rmid.*Mmarg_DDNeil{t}');
        sigma2_M_DDNeil(t)=(sum(M_Rmid.^2.*Mmarg_DDNeil{t}')-mu_M_DDNeil(t)^2);
        Iatc=find(at_R{t}==atc);
        FailProb_DDNeil(t)=sum(atmarg_DDNeil{t}(Iatc:end));
        Beta_DDNeil(t)=-norminv(FailProb_DDNeil(t));      
    end
    if length(a0_R)-1>=Lat_Targ && length(at_R{1})-1>=Lat_Targ && length(da_R{1})-1>=Lda_Targ && length(M_R{1})-1>=LM_Targ
       break; 
    end
    Error_a0marg_DDNeil=KLError(a0_R,a0marg_DDNeil);
    SumError_a0marg_DDNeil=sum(Error_a0marg_DDNeil);
    ErrorSort_a0marg_DDNeil=sort(Error_a0marg_DDNeil,2,'descend');
    if NumOfSplit>length(a0_R)-3 && Lat_Targ-(length(a0_R)-1)>0
        NumOfSplit_a0=1;
    elseif Lat_Targ-(length(a0_R)-1)<=NumOfSplit;
        NumOfSplit_a0=Lat_Targ-(length(a0_R)-1);        
    else
        NumOfSplit_a0=NumOfSplit;
    end
    a0Add_Ind=ismember(Error_a0marg_DDNeil,ErrorSort_a0marg_DDNeil(1:NumOfSplit_a0));
    a0Add=find(a0Add_Ind==1);
%     if length(a0_R)-1>=Lat_Targ
%         NumOfMerge_a0=length(a0_R)-1-Lat_Targ+NumOfMerge;
%         a0Merg_Ind=ismember(Error_a0marg_DDNeil,ErrorSort_a0marg_DDNeil(end-NumOfMerge_a0+1:end));
%         a0Merg=find(a0Merg_Ind==1);
%     end
    
    for t=1:T
        Error_damarg_DDNeil{t}=KLError(da_R{t},damarg_DDNeil{t});
        Error_atmarg_DDNeil{t}=KLError(at_R{t},atmarg_DDNeil{t});
        Error_Mmarg_DDNeil{t}=KLError(M_R{t},Mmarg_DDNeil{t});
        SumError_damarg_DDNeil(t)=sum(Error_damarg_DDNeil{t});
        SumError_atmarg_DDNeil(t)=sum(Error_atmarg_DDNeil{t});
        SumError_Mmarg_DDNeil(t)=sum(Error_Mmarg_DDNeil{t});
        ErrorSort_damarg_DDNeil=sort(Error_damarg_DDNeil{t},2,'descend');
        if NumOfSplit>length(da_R{t})-3 && Lda_Targ-(length(da_R{t})-1)>0
            NumOfSplit_da=1;
        elseif Lda_Targ-(length(da_R{t})-1)<=NumOfSplit;
            NumOfSplit_da=Lda_Targ-(length(da_R{t})-1);            
        else
            NumOfSplit_da=NumOfSplit;
        end      
        daAdd_Ind=ismember(Error_damarg_DDNeil{t},ErrorSort_damarg_DDNeil(1:NumOfSplit_da));
        daAdd{t}=find(daAdd_Ind==1);
%         if length(da_R{t})-1>=Lda_Targ
%             NumOfMerge_da=length(da_R{t})-1-Lda_Targ+NumOfMerge;
%             daMerg_Ind=ismember(Error_damarg_DDNeil{t},ErrorSort_damarg_DDNeil(end-NumOfMerge_da+1:end));
%             daMerg{t}=find(daMerg_Ind==1);
%         end
        ErrorSort_atmarg_DDNeil=sort(Error_atmarg_DDNeil{t},2,'descend');
        if NumOfSplit>length(at_R{t})-3 && Lat_Targ-(length(at_R{t})-1)>0
            NumOfSplit_at=1;
        elseif Lat_Targ-(length(at_R{t})-1)<=NumOfSplit;
            NumOfSplit_at=Lat_Targ-(length(at_R{t})-1);
        else
            NumOfSplit_at=NumOfSplit;
        end        
        atAdd_Ind=ismember(Error_atmarg_DDNeil{t},ErrorSort_atmarg_DDNeil(1:NumOfSplit_at));
        atAdd{t}=find(atAdd_Ind==1);
%         if length(at_R{t})-1>=Lat_Targ
%             NumOfMerge_at=length(at_R{t})-1-Lat_Targ+NumOfMerge;
%             atMerg_Ind=ismember(Error_atmarg_DDNeil{t},ErrorSort_atmarg_DDNeil(end-NumOfMerge_at+1:end));
%             atMerg{t}=find(atMerg_Ind==1);
%         end
        ErrorSort_Mmarg_DDNeil=sort(Error_Mmarg_DDNeil{t},2,'descend');
        if NumOfSplit>length(M_R{t})-3 && LM_Targ-(length(M_R{t})-1)>0
            NumOfSplit_M=1;
        elseif LM_Targ-(length(M_R{t})-1)<=NumOfSplit;
            NumOfSplit_M=LM_Targ-(length(M_R{t})-1);           
        else
            NumOfSplit_M=NumOfSplit;
        end                
        MAdd_Ind=ismember(Error_Mmarg_DDNeil{t},ErrorSort_Mmarg_DDNeil(1:NumOfSplit_M));
        MAdd{t}=find(MAdd_Ind==1); 
%         if length(M_R{t})-1>=LM_Targ
%             NumOfMerge_M=length(M_R{t})-1-LM_Targ+NumOfMerge;
%             MMerg_Ind=ismember(Error_Mmarg_DDNeil{t},ErrorSort_Mmarg_DDNeil(end-NumOfMerge_M+1:end)); 
%             MMerg{t}=find(MMerg_Ind==1);
%         end
        % Update the discretization according to saved add and merge
        % interval
        da_R{t}=DDStep(daAdd{t}, daMerg{t},da_R{t});
        at_R{t}=DDStep(atAdd{t}, atMerg{t},at_R{t});        
        M_R{t}=DDStep(MAdd{t}, MMerg{t},M_R{t});             
    end
    a0_R=DDStep(a0Add,a0Merg,a0_R);
    

% 
%     figure(1)
%     subplot(TotalIter,5,Iter*5-4);
%     DDCheck(da_R{1});
%     subplot(TotalIter,5,Iter*5-3);
%     DDCheck(da_R{5});    
%     subplot(TotalIter,5,Iter*5-2);
%     DDCheck(da_R{10}); 
%     subplot(TotalIter,5,Iter*5-1);
%     DDCheck(da_R{15});  
%     subplot(TotalIter,5,Iter*5);
%     DDCheck(da_R{20}); 
%     figure(2)
%     subplot(TotalIter,5,Iter*5-4);
%     DDCheck(at_R{1});
%     subplot(TotalIter,5,Iter*5-3);
%     DDCheck(at_R{5});    
%     subplot(TotalIter,5,Iter*5-2);
%     DDCheck(at_R{10}); 
%     subplot(TotalIter,5,Iter*5-1);
%     DDCheck(at_R{15});  
%     subplot(TotalIter,5,Iter*5);
%     DDCheck(at_R{20}); 
%     figure(3)
%     subplot(TotalIter,5,Iter*5-4);
%     DDCheck(M_R{1});
%     subplot(TotalIter,5,Iter*5-3);
%     DDCheck(M_R{5});    
%     subplot(TotalIter,5,Iter*5-2);
%     DDCheck(M_R{10}); 
%     subplot(TotalIter,5,Iter*5-1);
%     DDCheck(M_R{15});  
%     subplot(TotalIter,5,Iter*5);
%     DDCheck(M_R{20});    
%     figure(4)
%     subplot(TotalIter,1,Iter);
%     DDCheck(a0_R);

end

save('FailProb_DDNeil','FailProb_DDNeil','mu_at_DDNeil','sigma2_at_DDNeil','mu_da_DDNeil','sigma2_da_DDNeil');















