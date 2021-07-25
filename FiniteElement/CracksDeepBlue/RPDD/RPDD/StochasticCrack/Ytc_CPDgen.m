function Ytc_CPD=Ytc_CPDgen(alphaY_R,Y_R,Y_Rinf,Y_mean,Y_std,Y_mu,Y_sigma,DeltaN)
%disp('%%%%%%%%%%%%%%%%%% Ytc node CPD generation %%%%%%%%%%%%%%%%%%%');
LY=length(Y_R)-1;
LalphaY=length(alphaY_R)-1;
Ytc_CPD=zeros(LY,LalphaY,LY);
for i=1:LalphaY
    for j=1:LY
        Ytc_CDF=zeros(1,LY-1);
        for k=1:LY-1
            Count=0;
            Ytc_cdf=0;
            if j~=LY
                for ialphaY=(alphaY_R(i+1)+alphaY_R(i))/2;
                %for ialphaY=[alphaY_R(i) (alphaY_R(i)+alphaY_R(i+1))/2 alphaY_R(i+1)]
                    for jYtm=(Y_R(j)+Y_R(j+1))/2;
                    %for jYtm=[Y_R(j) (Y_R(j)+Y_R(j+1))/2 Y_R(j+1)]
                        Ytc_cdf=Ytc_cdf+fYint(ialphaY,jYtm,Y_R(k+1),Y_mean,Y_std,Y_mu,Y_sigma,DeltaN);
                        Count=Count+1;
                    end
                end
                         
            else
                for ialphaY=(alphaY_R(i+1)+alphaY_R(i))/2;
                %for ialphaY=[alphaY_R(i) (alphaY_R(i)+alphaY_R(i+1))/2 alphaY_R(i+1)]
                    for jYtm=Y_R(j);
                    %for jYtm=[Y_R(j) 1000]
                        Ytc_cdf=Ytc_cdf+fYint(ialphaY,jYtm,Y_R(k+1),Y_mean,Y_std,Y_mu,Y_sigma,DeltaN);
                        Count=Count+1;
                    end
                end
            end
            Ytc_cdf=Ytc_cdf/Count;
            Ytc_CDF(k)=Ytc_cdf;
        end 
        Ytc_CDF2=[0 Ytc_CDF(1:end-1)];
        Ytc_pdf=Ytc_CDF-Ytc_CDF2;
        if sum(Ytc_pdf<0)>=1
           disp('warning: pdf of Y is negative') 
           [i,j,find(Ytc_pdf<0)]
           Ytc_CDF
           Ytc_pdf(find(Ytc_pdf<0))
           Ytc_pdf(find(Ytc_pdf<0))=zeros(1,length(find(Ytc_pdf<0)));
           Ytc_pdf
        end
        if sum(Ytc_pdf>1)>=1
            disp('warning:pdf of Y is larger than one')
            [i,j,find(Ytc_pdf>1)]
            Ytc_pdf(find(Ytc_pdf<0))
        end
        if sum(Ytc_pdf)<0 || sum(Ytc_pdf)>1
            disp('warning: the prob. of last region of Ytm is negative')
            sum(Ytc_pdf)
            Ytc_pdf=[Ytc_pdf 0];
        else
            Ytc_pdf=[Ytc_pdf (1-sum(Ytc_pdf))];
        end
        Ytc_CPD(j,i,:)=Ytc_pdf;
    end
end