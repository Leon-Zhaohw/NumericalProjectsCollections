function [damarg_myDD, atmarg_myDD, Mmarg_myDD, da_R, at_R, M_R]=...
    ConstCrack_myDDKL (mu_a0, sigma_a0, mu_da, sigma_da, sigma_error,Evidence, atc, T,da_R,at_R,M_R,Lda_Targ, Lat_Targ)

%Constant crack growth model with RPDD algorithm and KL error judgement. 


TotalIter=20;

for Iter=1:TotalIter
    % Intialize the Add Merge Interval for at and da

    Lda=length(da_R)-1;
    Lat=length(at_R)-1;
    LM=length(M_R)-1;


    [a0_CPD,a1_CPD, da_CPD, at_CPD,M_CPD]=CPDgen_DBN(at_R,mu_a0,sigma_a0,da_R,mu_da, M_R, sigma_da,sigma_error,atc);

    
    ev=cell(1,T*3);
    SizeEV=size(Evidence);
    for i=1:SizeEV(1)
        for j=1:LM
            if M_R(j+1)>=Evidence(i,2)
                ev{1,Evidence(i,1)*3}=j;
                break
            end
        end
    end  
    % send to inference engin
    [SBNengine,damarg_myDD, atmarg_myDD, Mmarg_myDD]=ConstCrack_DBNModel(T,ev,da_CPD,a1_CPD,at_CPD,M_CPD);

    
    if Lda>=Lda_Targ && Lat>=Lat_Targ
        break
    end
    
    %%%%%%%%%%%%%%%%%%%%% Stoping Judgment %%%%%%%%%%%%%%%%%%%%%%%    
    %     if Iter==TotalIter_myDD
    %         break
    %     end
    %     if Iter>1
    %         StopJudge(Iter)=sum(abs(Beta_myDD(:,Iter)-Beta_myDD(:,Iter-1)))/sum(abs(Beta_myDD(:,Iter)))/T
    %         StopJudge_Smooth(Iter)=sum(abs(Beta_myDD(1:max(Evidence(:,1)),Iter)...
    %             -Beta_myDD(1:max(Evidence(:,1)),Iter-1)))/sum(abs(Beta_myDD(1:max(Evidence(:,1)),Iter)))/length(1:max(Evidence(:,1)))
    %         StopJudge_Predict(Iter)=sum(abs(Beta_myDD(max(Evidence(:,1))+1:end,Iter)-...
    %             Beta_myDD(max(Evidence(:,1))+1:T,Iter-1)))/sum(abs(Beta_myDD(max(Evidence(:,1))+1:T,Iter)))/length(max(Evidence(:,1))+1:T)
    %         if StopJudge(Iter)<=0.005
    %             break
    %         end
    %     end
    %%%%%%%%%%% Generate the posterior distribution after inserting virtual evidence
     Iatc=find(at_R==atc);
     Count1=1;
     damarg2_myDD=cell(1);
     atmarg2_myDD=cell(1);
     for t=1:max(Evidence(:,1))
         Count2=1;
         Count3=1;
         ev2=ev;
         % Insert virtual evidence at first non-zeros bin of posterior
         % distribution 
         for VirtualEv=Lat:-1:1
             if atmarg_myDD(VirtualEv,t)>0
                 break
             end
         end
         ev2{1,t*3-1}=VirtualEv;
         [engine_ev,loglik_ev]=enter_evidence(SBNengine,ev2);
         atmarg2_myDD{Count1}=zeros(Lat,1);
         damarg2_myDD{Count1}=zeros(Lda,1);
         for tt=1:t-1
             marg_ev=marginal_nodes(engine_ev,tt*3-1);
             atmarg2_myDD{Count1}(:,Count2)=marg_ev.T;
             Count2=Count2+1;
         end
         for tt=1:max(Evidence(:,1))
             marg_ev=marginal_nodes(engine_ev,tt*3-2);
             damarg2_myDD{Count1}(:,Count3)=marg_ev.T;  
             I=find(max(damarg2_myDD{Count1}(:,Count3))==damarg2_myDD{Count1}(:,Count3));
             Count3=Count3+1;
         end
         Count1=Count1+1;
     end

    %%%%%%%%%%%%% Generate the add and merge interval for da and at nodes %%%%%%%%%%%%%%
    atMerg=[];
    daMerg=[];
    MMerg=[];
    atAdd=[];
    daAdd=[];
    MAdd=[];
    if 1
        %disp('%%%%%%%%%%%%%%%%%%%%%%%% Add Step %%%%%%%%%%%%%%%%%%%%%%%%%')
        daAdd1=AddStep_KL(damarg_myDD,da_R);
        daAdd2=AddStep_KL(damarg2_myDD,da_R);
        atAdd1=AddStep_KL(atmarg_myDD,at_R);
        atAdd2=AddStep_KL(atmarg2_myDD,at_R); 

        atAdd=[atAdd1 atAdd2];
        atAdd=sort(unique(atAdd),'ascend');
        % the last interval is betwen atc to inf which will nerver give
        % subdividion on that bin
        if atAdd(end)==Lat
            atAdd(end)=[];
        end
        if length(atAdd)+length(at_R)-1>=Lat_Targ
            atAdd=randsample(atAdd,Lat_Targ-length(at_R)+1);
        end     
        
        daAdd=[daAdd1 daAdd2];
        daAdd=sort(unique(daAdd),'ascend');
        % the last interval is betwen daLimit to inf which will nerver give
        % subdividion on that interval
        if daAdd(end)==Lda
            daAdd(end)=[];
        end 
        if length(daAdd)+length(da_R)-1>=Lda_Targ
            daAdd=randsample(daAdd,Lda_Targ-length(da_R)+1);
        end       
        if isempty(atAdd)
            MAdd=[];
        else
            MAdd=atAdd;
            for i=length(MAdd):-1:1
                if MAdd(i)>Iatc
                    MAdd(i)=[];
                end            
            end
            if MAdd(end)==LM
                MAdd(end)=[];
            end 
        end

    end
    % Update the discretization according to saved add and merge interval  
    
    da_R=DDStep(daAdd, daMerg,da_R);
    at_R=DDStep(atAdd, atMerg,at_R);
    M_R=DDStep(MAdd, MMerg,M_R);  

end














