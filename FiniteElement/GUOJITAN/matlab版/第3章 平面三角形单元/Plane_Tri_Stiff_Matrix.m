
function ZK= Plane_Tri_Stiff_Matrix
% 计算结构总刚度矩阵
%调用函数： 弹性矩阵函数Plane_Elastic_Matrix,  几何矩阵函数Plane_B_Matrix。
%计算单元刚度矩阵并组装为总刚度矩阵
%返回量：总刚度矩阵ZK(2*nd,2*nd )
global  pm  E  nv  t0  nd  ne  XY  EL
%全局变量：问题性质,弹性模量,泊松比,厚度、结点数、单元数、结点坐标、单元信息
ZK = zeros(2*nd,2*nd ) ;            % 结构的总刚度矩阵
D = Elastic_Matrix (pm, E, nv);        %调用函数的弹性矩阵，E-弹性模量，nv—泊松比
for ie=1:1:ne                                     %对单元循环
    [B,A] = Plane_B3_Matrix( ie) ;                    %调用函数，计算单元ie的几何矩阵及面积
    S = D * B ;                                    %应力矩阵
    KE = t0*A*transpose(S)*B;                         %单元刚度矩阵
    %  将单元刚度矩阵KE集成到整体刚度矩阵ZK
    for r=1:1:3
        i0=2* EL(ie,r);
        m0 = 2*r ;
        for s=1:1:3
            j0=2* EL(ie,s);
            n0 =2*s;
            %将单刚中与r、s相对应的2×2子阵，叠加到总刚中
            ZK([i0-1,i0],[j0-1,j0]) = ZK([i0-1,i0],[j0-1,j0]) + KE([m0-1,m0],[n0-1,n0]) ;
        end
    end
end
return






