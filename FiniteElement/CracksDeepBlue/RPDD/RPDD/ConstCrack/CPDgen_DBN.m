function [a0_CPD,a1_CPD,da_CPD, at_CPD,M_CPD]=CPDgen(at_R,mu_a0,sigma_a0,da_R,mu_da, M_R, sigma_da,sigma_error,atc)


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
a1_CPD=zeros(Lda,Lat);
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
a_Rinf=100;
M_CPD=zeros(Lat, LM);
for i=1:Lat
   M_PDF=zeros(1,LM-1);
   if i~=Lat
       for j=1:LM-1
           if j==1
               M_PDF(j)=normcdf(M_R(j+1),(at_R(i)+at_R(i+1))/2,sigma_error);
           else
               M_PDF(j)=normcdf(M_R(j+1),(at_R(i)+at_R(i+1))/2,sigma_error)-...
                        normcdf(M_R(j),(at_R(i)+at_R(i+1))/2,sigma_error);
           end
       end
   else
       for j=1:LM-1
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
   if i~=Lat && M_PDF(end-1)==0 && M_PDF(end)~=0
       disp('warning: the pdf is not smooth')
       i
       log10(M_PDF(end))
       M_PDF(end)=0;
   end
    M_CPD(i,:)=M_PDF;
end  