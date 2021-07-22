

function  Plane_Triangular_Element
% 本程序为采用线性三角形单元求解平面问题，计算在自重等体力、集中力以及
%线性分布面力作用下的变形和应力，并将结果存储到文件。
%计算结果存储在： sjx_RES.txt
% 调用5个功能函数：有限元模型数据、计算结构总刚、载荷向量、求解有限元方程、计算应力
%   [file_in,file_out] = File_Name      %输入文件名及计算结果输出文件名
file_out='sjx_RES.txt'
Plane_Tri_Model_Data  ;                 %有限元模型数据
%   Plane_Tri_Modle_Figure;             %显示有限元模型图形，以便于检查
ZK= Plane_Tri_Stiff_Matrix              %计算结构总刚
ZQ= Plane_Tri_Load_Vector                %计算总的载荷向量
U= Solve_1_Model(file_out,2,ZK,ZQ);     %求解有限元方程，得到结点位移，保存到文件并返出
Stress_nd=Plane_Tri_Strees(file_out,U);   %应力分析，将计算结果保存到文件，并返出绕点平均应力
%    Plane_Tri_Post_Contour(U,Stress_nd );     %后处理模块，显示变形图、不同应力分量的云图
fclose all;
end
