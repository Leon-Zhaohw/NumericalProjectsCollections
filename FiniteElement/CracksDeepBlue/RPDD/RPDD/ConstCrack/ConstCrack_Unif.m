% This is a comparision plot for constant crack with uniform discretization
% This code should be run after the "ConsCrack_myDD.m"

function [damarg_Unif, atmarg_Unif, Mmarg_Unif]=ConstCrack_Unif(mu_a0, sigma_a0, mu_da, sigma_da, sigma_error,Evidence, atc, T,da_R, at_R,M_R)

Lda=length(da_R)-1;
Lat=length(at_R)-1;
LM=length(M_R)-1;
%CPT table of a0
a0_cdf=normcdf(at_R,mu_a0,sigma_a0);
a0_cdf2=a0_cdf(2:end);
a0_CPD=a0_cdf2-a0_cdf(1:end-1);
%CPT table of da
da_cdf=normcdf(da_R,mu_da,sigma_da);
da_cdf2=da_cdf(2:end);
da_CPD=da_cdf2-da_cdf(1:end-1);
%CPT table of at
at_CPD=zeros(Lat,Lda,Lat);
for i=1:Lat
    for j=1:Lda
        Count=0;
        for iatm=[at_R(i) (at_R(i+1)+at_R(i))/2 at_R(i+1)]
            for jda=[da_R(j) (da_R(j+1)+da_R(j))/2 da_R(j+1)]
                at_approx=iatm+jda;
                for k=1:Lat 
                    if at_R(k+1)>=at_approx
                        at_CPD(i,j,k)=at_CPD(i,j,k)+1;
                        break
                    end
                end
                Count=Count+1;
            end
        end
        at_CPD(i,j,:)=at_CPD(i,j,:)/Count;
    end
end
%CPT table of a1
a1_CPD=zeros(Lda ,Lat );
for i=1:Lda 
    for j=1:Lat 
        Count=0;
        a0_prob=a0_CPD(j);
        for ida=[da_R(i) (da_R(i)+da_R(i+1))/2 da_R(i+1)]
            for ja0=[at_R(j) (at_R(j)+at_R(j+1))/2 at_R(j+1)]
                at_approx=ja0+ida;
                for k=1:Lat 
                   if at_R(k+1)>=at_approx
                       a1_CPD(i,k)=a1_CPD(i,k)+1*a0_prob;
                       break
                   end
                end
                Count=Count+1;
            end
        end
        a1_CPD(i,j,:)=a1_CPD(i,j,:)/Count;
    end
end
%CPT table of M node
lambda=1;
a_Rinf=1000;
M_CPD=zeros(Lat , LM );
for i=1:Lat 
   M_PDF=zeros(1,LM -1);
   if i~=Lat 
       for j=1:LM -1
           if j==1
               M_PDF(j)=normcdf(M_R(j+1),(at_R(i)+at_R(i+1))/2,sigma_error);
           else
               M_PDF(j)=normcdf(M_R(j+1),(at_R(i)+at_R(i+1))/2,sigma_error)-...
                        normcdf(M_R(j),(at_R(i)+at_R(i+1))/2,sigma_error);
           end
       end
   else
       for j=1:LM -1
           if j==1
               M_PDF(j)=quadgk(@(at)normcdf(M_R(j+1),at,sigma_error).*lambda.*exp(-lambda.*(at-atc)),at_R(i),a_Rinf);
           else
               M_PDF(j)=quadgk(@(at)normcdf(M_R(j+1),at,sigma_error).*lambda.*exp(-lambda.*(at-atc)),at_R(i),a_Rinf)-...
                        quadgk(@(at)normcdf(M_R(j),at,sigma_error).*lambda.*exp(-lambda.*(at-atc)),at_R(i),a_Rinf);
           end
       end
   end
   % CPT check
   if sum(M_PDF)<0 || (sum(M_PDF)-1)>0
        disp('warning: the pdf of M may be not right')
        i
        1-sum(M_PDF)
        M_PDF=[M_PDF 0];
   else
        M_PDF=[M_PDF 1-sum(M_PDF)];
   end
   if i~=LM  && M_PDF(end-1)==0 && M_PDF(end)~=0
       disp('warning: the pdf is not smooth')
       i
       log10(M_PDF(end))
       M_PDF(end)=0;
   end
    M_CPD(i,:)=M_PDF;
end           


%Build the Bayesian Network
NN=T*3;
hnodes1=1:3:3*T;
hnodes2=2:3:3*T;
onodes=3:3:3*T;

eclass=ones(1,NN);

H1class=1;
H2class=2;
H3class=4;
O1class=3;

eclass(hnodes1)=H1class;
eclass(hnodes2(1))=H2class;
eclass(hnodes2(2:end))=H3class;
eclass(onodes(1:end))=O1class;
% Define the linkage relation
dag=zeros(NN);
for i=1:T
    dag(i*3-1,i*3)=1;
    dag(i*3-2,i*3-1)=1;
end
for i=1:T-1
    dag(i*3-1,i*3+2)=1;
end
% Define the nodal size
ns=zeros(NN,1);
for i=1:3:NN
    ns(i)=Lda ;
    ns(i+1)=Lat ;
    ns(i+2)=LM ;
end

dnodes=1:NN;
bnet=mk_bnet(dag,ns,'discrete',dnodes,'equiv_class',eclass);

bnet.CPD{H1class}=tabular_CPD(bnet,hnodes1(1),da_CPD);

bnet.CPD{H2class}=tabular_CPD(bnet,hnodes2(1),a1_CPD);

bnet.CPD{H3class}=tabular_CPD(bnet,hnodes2(2),at_CPD);

bnet.CPD{O1class}=tabular_CPD(bnet,onodes(1),M_CPD);
% Inference Engine
SBNengine = jtree_inf_engine(bnet);
% Generate the evidence
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
[engine_ev,loglik_ev]=enter_evidence(SBNengine,ev);
% Generate the posterior distribution and store it
damarg_Unif=zeros(Lda ,T);
atmarg_Unif=zeros(Lat ,T);
Mmarg_Unif=zeros(LM ,T);
for t=1:T
    marg_ev=marginal_nodes(engine_ev,t*3-2);
    damarg_Unif(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(engine_ev,t*3-1);
    atmarg_Unif(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(engine_ev,t*3);
    if length(marg_ev.T)==1
        marg_ev.T=zeros(LM ,1);
    end
    Mmarg_Unif(:,t)=marg_ev.T;
end


% figure
% plot(Lat,log10(BetaError_DD),'r','linewidth',3);
% hold on
% plot(Lat,log10(BetaError_Unif),'--b','linewidth',3);
% legend('RPDD','uniform');
% xlabel('number of intervals for a_t','fontsize',15)
% ylabel('log_1_0(\epsilon_\beta)','fontsize',15);
% set(gca,'fontsize',15)
% grid on 
% 
% 





