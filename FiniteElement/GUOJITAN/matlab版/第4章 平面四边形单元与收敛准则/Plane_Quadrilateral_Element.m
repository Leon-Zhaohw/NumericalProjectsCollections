

function  Plane_Quadrilateral_Element
% 本程序为采用四点等参元求解平面问题，计算在自重等体力、集中力、线性分布面力作用下的变形和应力，并将结果存储到文件。
%计算结果存储在： sddcy_RES.txt
% 调用5个功能函数：有限元模型数据、计算结构总刚、载荷向量、求解有限元方程、计算应力
%  [file_in,file_out] = File_Name ;  %输入文件名及计算结果输出文件名
file_out='Q4dengcanyuan_RES.txt'
Plane_Quad4_Model;              % 读入有限元模型数据并图形显示
%    Plane_Q4_Modle_Figure;                 %显示有限元模型图形，以便于检查
ZK= Plane_Quad_4_Stiff_Matrix;            %计算结构总刚
ZQ= Plane_Quad_4_Load_Matrix;          %计算总的载荷向量
U= Solve_1_Model(file_out,2,ZK,ZQ);     %求解有限元方程，得到结点位移，保存到文件并返出
Stress_nd = Quadrilateral_Strees(file_out, U);        %应力分析，并将计算结果保存到文件中
%  Plane_Quad_4_Post_Contour(U,Stress_nd );         %后处理模块，显示变形图、不同应力分量的云图
fclose all;
end












