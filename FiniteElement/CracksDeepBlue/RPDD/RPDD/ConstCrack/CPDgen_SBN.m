function  [a0_CPD,da_CPD, at_CPD,M_CPD]=CPDgen_SBN(a0_R,da_R,at_R,M_R,mu_a0,sigma_a0,mu_da, sigma_da,sigma_error,atc)


%CPT table of a0

a0_cdf=normcdf(a0_R,mu_a0,sigma_a0);
a0_cdf2=a0_cdf(2:end);
a0_CPD=a0_cdf2-a0_cdf(1:end-1);

T=length(da_R);
for t=1:T
    Lda=length(da_R{t})-1;
    Lat=length(at_R{t})-1;
    LM=length(M_R{t})-1;
    if t>1
        Latm=length(at_R{t-1})-1;
        atm_R=at_R{t-1};
    else
        Latm=length(a0_R)-1;
        atm_R=a0_R;
    end
    
    %CPT table of da
    da_cdf=normcdf(da_R{t},mu_da,sigma_da);
    da_cdf2=da_cdf(2:end);
    da_CPD{t}=da_cdf2-da_cdf(1:end-1);

    %CPT table of at
    at_CPD{t}=zeros(Latm,Lda,Lat);
    for i=1:Latm
        for j=1:Lda
            Count=0;
            for iatm=[atm_R(i):(atm_R(i+1)-atm_R(i))/2:atm_R(i+1)]
                for jda=[da_R{t}(j):(da_R{t}(j+1)-da_R{t}(j))/2:da_R{t}(j+1)]
                    at_approx=iatm+jda;
                    for k=1:Lat
                        if at_R{t}(k+1)>=at_approx
                            at_CPD{t}(i,j,k)=at_CPD{t}(i,j,k)+1;
                            break
                        end
                    end
                    Count=Count+1;
                end
            end
            at_CPD{t}(i,j,:)=at_CPD{t}(i,j,:)/Count;
        end
    end

    %CPT table of M node
    lambda=1;
    a_Rinf=100;
    M_CPD{t}=zeros(Lat, LM);
    for i=1:Lat
       M_PDF=zeros(1,LM-1);
       if i~=Lat
           for j=1:LM-1
               if j==1
                   M_PDF(j)=normcdf(M_R{t}(j+1),(at_R{t}(i)+at_R{t}(i+1))/2,sigma_error);
               else
                   M_PDF(j)=normcdf(M_R{t}(j+1),(at_R{t}(i)+at_R{t}(i+1))/2,sigma_error)-...
                            normcdf(M_R{t}(j),(at_R{t}(i)+at_R{t}(i+1))/2,sigma_error);
               end
           end
       else
           for j=1:LM-1
               if j==1
                   M_PDF(j)=quadgk(@(at)normcdf(M_R{t}(j+1),at,sigma_error).*lambda.*exp(-lambda.*(at-atc)),at_R{t}(i),a_Rinf);
               else
                   M_PDF(j)=quadgk(@(at)normcdf(M_R{t}(j+1),at,sigma_error).*lambda.*exp(-lambda.*(at-atc)),at_R{t}(i),a_Rinf)-...
                            quadgk(@(at)normcdf(M_R{t}(j),at,sigma_error).*lambda.*exp(-lambda.*(at-atc)),at_R{t}(i),a_Rinf);
               end
           end
       end
       % CPT check
       if sum(M_PDF)<0 || (sum(M_PDF)-1)>0
            disp('warning: the pdf of M may be not right')
            at_R{t}
            M_R{t}
            M_PDF
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
        M_CPD{t}(i,:)=M_PDF;
    end 
end