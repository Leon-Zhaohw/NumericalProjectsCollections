function [SBNengine alphaY_CPD Y1_CPD Ytc_CPD a1_CPD at_CPD M_CPD]=Crack_SBNengine_M(T,alphaY_R,Y_R,at_R,M_R,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc)

[alphaY_CPD Y1_CPD Ytc_CPD a1_CPD,at_CPD,M_CPD]=CPDgen(alphaY_R,Y_R,at_R,M_R,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc);

NN=T*4;       %Num of total nodes
hnodes1=1:4:NN;
hnodes2=2:4:NN;
hnodes3=3:4:NN;
onodes=4:4:NN;
H1class=1;
H2class=2;
H3class=3;
H4class=4;
H5class=5;
O1class=6;

eclass=ones(1,NN);
eclass(onodes)=O1class;   %M
eclass(hnodes1)=H1class;  %alphaY
eclass(hnodes2(1))=H2class;  %Y1
eclass(hnodes3(1))=H3class;  %a1
eclass(hnodes2(2:end))=H4class; %>=Y2
eclass(hnodes3(2:end))=H5class;  %>=a2

dag=zeros(NN);

for i=1:4:NN
   dag(i,i+1)=1;
   dag(i+1,i+2)=1;
   dag(i+2,i+3)=1;
end
for i=1:4:NN-4
   dag(i+1,i+5)=1;
   dag(i+2,i+6)=1;
end

ns=zeros(NN,1);

LalphaY=length(alphaY_R)-1;
LY=length(Y_R)-1;
Lat=length(at_R)-1;
LM=length(M_R)-1;
for i=1:4:NN
    ns(i)=LalphaY;
    ns(i+1)=LY;
    ns(i+2)=Lat;
    ns(i+3)=LM;
end


dnodes=1:NN;
bnet=mk_bnet(dag,ns,'discrete',dnodes,'equiv_class',eclass);

bnet.CPD{H1class}=tabular_CPD(bnet,hnodes1(1),alphaY_CPD);

bnet.CPD{H2class}=tabular_CPD(bnet,hnodes2(1),Y1_CPD);

bnet.CPD{H3class}=tabular_CPD(bnet,hnodes3(1),a1_CPD);

bnet.CPD{H4class}=tabular_CPD(bnet,hnodes2(2),Ytc_CPD);

bnet.CPD{H5class}=tabular_CPD(bnet,hnodes3(2),at_CPD);

bnet.CPD{O1class}=tabular_CPD(bnet,onodes(1),M_CPD);

SBNengine = jtree_inf_engine(bnet);


