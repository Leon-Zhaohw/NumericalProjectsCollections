

function  ZK= Plane_Quad_4_Stiff_Matrix
% 高斯积分法计算四点等参元的结构总刚度矩阵
%需要调用的函数：（1）弹性矩阵D， （2）四点等参元单刚矩阵。
%输入量：结点坐标、单元信息，用全局变量赋值
%返回量：总刚度矩阵ZK(2*nd,2*nd )
global  pm  E  nv  t0  nd  ne  XY  EL
%全局变量：问题性质,弹性模量,泊松比,厚度、结点数、单元数、结点坐标、单元信息、弹性矩阵
ZK = zeros(2*nd,2*nd ) ;            % 结构的总刚度矩阵
D =Elastic_Matrix (pm, E, nv);
for ie=1:1:ne                                     %对单元循环
    KE =Gauss_Stiff_Matrix( ie,D );                        %单元刚度矩阵
    %  将单元刚度矩阵KE集成到整体刚度矩阵ZK
    for r=1:1:4
        i0=2* EL(ie,r);
        m0 = 2*r;
        for s=1:1:4
            j0=2* EL(ie,s);
            n0 =2*s;
            %将单刚中与r、s相对应的2×2子阵，叠加到总刚中
            ZK((i0-1):i0, (j0-1):j0) = ZK((i0-1):i0, (j0-1):j0)+ KE((m0-1):m0, (n0-1):n0) ;
        end
    end
end
return














