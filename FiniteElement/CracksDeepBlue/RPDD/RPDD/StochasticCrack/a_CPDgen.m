function [at_CPD,a1_CPD]=a_CPDgen(a0_CPD,at_R,Y_R,C,m,DeltaS,DeltaN)
disp('%%%%%%%%%%%%%%%%%%% at node CPD generation %%%%%%%%%%%%%%%%%%%%%%%')
Lat=length(at_R)-1;
LY=length(Y_R)-1;
at_CPD=zeros(Lat,LY,Lat);
a1_CPD=zeros(LY,Lat);
for i=1:LY
    for j=1:Lat
        a1_tmp=zeros(Lat,1);
        at_tmp=zeros(Lat,1);
        count=0;
        if j~=Lat
            if i~=LY
                for jat=[at_R(j) (at_R(j+1)+at_R(j))/2  at_R(j+1)]
                %for jat=at_R(j):(at_R(j+1)-at_R(j))/50:at_R(j+1)
                    for iY=(Y_R(i)+Y_R(i+1))/2;
                    %for iY=(Y_R(i)+0.01):(Y_R(i+1)-Y_R(i)-0.02)/50:Y_R(i+1)-0.01
                        a_act=faint(jat,iY,C,m,DeltaS,DeltaN);
                        count=count+1; 
                        for k=1:Lat
                            if at_R(k+1)>=a_act
                                at_tmp(k)=at_tmp(k)+1;
                                a1_tmp(k)=a1_tmp(k)+a0_CPD(j)*1; 
                                break
                            end
                        end  
                    end
                end
            else
                for jat=[at_R(j)  (at_R(j+1)+at_R(j))/2  at_R(j+1)] 
                %for jat=at_R(j):(at_R(j+1)-at_R(j))/100:at_R(j+1)
                    for iY=Y_R(i);  
                    %for iY=[Y_R(i) 100]
                        a_act=faint(jat,iY,C,m,DeltaS,DeltaN);
                        count=count+1;
                        for k=1:Lat
                            if at_R(k+1)>=a_act
                                at_tmp(k)=at_tmp(k)+1;
                                a1_tmp(k)=a1_tmp(k)+a0_CPD(j)*1;
                                break
                            end
                        end
                    end
                end
            end
            at_tmp=at_tmp/count;
            a1_tmp=a1_tmp/count;
            if a_act<=at_R(j)
               disp('warning: the at is less than atm')
                %[i,j,k]
                %a_act
            end
        else
            at_tmp(Lat)=1;
            a1_tmp(Lat)=a1_tmp(Lat)+1*a0_CPD(Lat);
        end 
        
        if abs(sum(at_tmp)-1)>=10^(-8)
             i,j
             sum(at_tmp)
           disp('warning of at_tmp') 
        end
        a1_CPD(i,:)=a1_CPD(i,:)+a1_tmp';
        at_CPD(j,i,:)=at_tmp;    
    end
end

for i=1:LY
    if abs(sum(a1_CPD(i,:))-1)>=1E-8
        disp('warning probability of a1 maybe wrong')
    end
    for j=1:Lat
        if abs(sum(at_CPD(j,i,:))-1)>=1E-8
            disp('warning probability of at maybe wrong')
        end
    end
end

