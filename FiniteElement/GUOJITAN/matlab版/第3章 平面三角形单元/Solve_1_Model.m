

function U= Solve_1_Model(file_out,ndim,ZK,ZQ)
%  采用置“1”法解有限元方程，得到总的结点位移向量
% 输入参数：file_out,——与模型文件相应的txt计算结果输出文件名
% ndim-弹性体维数: 平面问题—2, 轴对称问题—2 , 空间问题—3；
% ZK----总刚， ZQ----总载荷向量
%  返回值：总的结点位移向量
global  pm  nd ne ng  BC   %全局变量pm—问题属性，nd—结点总数，ng—约束总数，BC数组—边界条件
U=zeros(ndim*nd,1);
[n1 n2]=size(BC);
for r=1:1:ng
    i0= ndim *(BC(r,1)-1) + BC(r,2);
    ZK(i0,:) =  0 ;  ZK(:,i0) =  0 ;
    ZK(i0,i0) =  1.0 ;
    ZQ(i0,1)=  0;
    if  n2==3   &  BC(r,3)~=0      % 给定的位移为非零值
        ZQ(i0,1)= BC(r,3) ;
    end
end
disp('正在求解有限元方程，请稍候......')
U = ZK\ZQ ;
ty={'平面应力问题','平面应变问题','轴对称问题','空间问题'};
fid=fopen(file_out,'wt');
fprintf(fid,'         有限元计算结果 \n');
fprintf(fid,' 问题属性：%s,   结点总数：%i，   单元总数：%i \n',ty{pm},nd,ne);
fprintf(fid,'    结点位移计算结果 \n') ;               % 设定输出到文件
fprintf(fid,'\n    Node        Ux          Uy     \n') ;
for j=1:1:nd                                  %将结点位移计算结果存储到文件
    fprintf(fid,'\n %8i%20.6e%20.6e\n',j, U(2*j-1,1),U(2*j,1));       
end
return











