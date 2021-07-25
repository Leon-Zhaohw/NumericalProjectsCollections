function M_CPD=M_CPDgen(at_R,M_R,M_sigma,lambda,ac,a_Rinf)
%disp('%%%%%%%%%%%%%%%%%%%%% M node CPD genration %%%%%%%%%%%%%%%%%%%%%%%%')
Lat=length(at_R)-1;
LM=length(M_R)-1;
M_CPD=zeros(Lat,LM);
for i=1:Lat
    M_PDF=zeros(1,LM-1);
    if i~=Lat
        for j=1:LM-1
            Count=0;
            %for iat=at_R(j):(at_R(j+1)-at_R(j))/100:at_R(j+1)
            for iat=[at_R(i) (at_R(i)+at_R(i+1))/2 at_R(i+1)]
                if j==1
                    M_PDF(j)=M_PDF(j)+normcdf(M_R(j+1),iat,M_sigma);
                else
                    M_PDF(j)=M_PDF(j)+normcdf(M_R(j+1),iat,M_sigma)-...
                            normcdf(M_R(j),iat,M_sigma);
                end
                Count=Count+1;
            end
            M_PDF(j)=M_PDF(j)/Count;
        end
    else
        for j=1:LM-1
            if j==1
                M_PDF(j)=quadgk(@(at)normcdf(M_R(j+1),at,M_sigma).*lambda.*exp(-lambda.*(at-ac)),at_R(i),a_Rinf);
            else
                M_PDF(j)=quadgk(@(at)normcdf(M_R(j+1),at,M_sigma).*lambda.*exp(-lambda.*(at-ac)),at_R(i),a_Rinf)-...
                quadgk(@(at)normcdf(M_R(j),at,M_sigma).*lambda.*exp(-lambda.*(at-ac)),at_R(i),a_Rinf)+...
                quadgk(@(at)normcdf(-M_R(j),at,M_sigma).*lambda.*exp(-lambda.*(at-ac)),at_R(i),a_Rinf)-...
                quadgk(@(at)normcdf(-M_R(j+1),at,M_sigma).*lambda.*exp(-lambda.*(at-ac)),at_R(i),a_Rinf);
            end
        end
    end
    if sum(M_PDF)<=0 || (sum(M_PDF)-1)>0
        disp('warning: the pdf of M may be not right')
        sum(M_PDF)
        M_PDF=[M_PDF 0];
        M_PDF=M_PDF/sum(M_PDF);
    else
        M_PDF=[M_PDF 1-sum(M_PDF)];
    end
    M_CPD(i,:)=M_PDF;
end



