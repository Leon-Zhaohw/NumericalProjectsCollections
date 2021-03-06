

function   Stress_nd = Quadrilateral_Strees(file_out, U)
%  求单元高斯点或结点应力分量及主应力、Mises应力
%  输入参数：U----总的结点位移向量； file_out—字符串，存储计算结果的文件名
%  返回值：无%
global  pm E  nv  t0  lou nd ne  XY  EL
%全局变量：结点数、单元数、结点坐标、单元信息 、弹性矩阵
fid=fopen(file_out,'at');      %打开存储计算结果的文件
%计算单元应力
St_nd=zeros(ne,12);
%选择计算应力的位置点：高斯点？ 结点？
%g_node=input('请选择计算应力的位置点：高斯积分点应力,请输（1）；结点应力,请输（2）；——') ;
for g_node=1:2
    if g_node ==1
        fprintf( '\n     正在计算单元高斯积分点应力...\n')
        fprintf(fid,'\n          单元高斯积分点应力分量计算结果 \n')          %高斯点表头
        fprintf(fid,'\n       1#（-0.57735，-0.57735）        2#（0.57735，-0.57735）         3#（0.57735， 0.57735）          4#（-0.57735，0.57735）\n')
        fprintf(fid,'  单元号    sigx  sigy   tau    sigx  sigy   tau    sigx   sigy   tau  sigx   sigy    tau \n')
        gx = [-0.577350269189626,  0.577350269189626] ;          %高斯点坐标
    else
        fprintf( '\n     正在计算单元结点处应力...\n')
        fprintf(fid,'\n          四点等参元结点应力分量\n')                  %结点表头
        fprintf(fid,'\n       1#（-1，-1）      2#（1，-1）       3#（1， 1）  4#（-1，1 ）\n')
        gx = [-1,  1] ;                                         %结点坐标
        fprintf(fid,'  结点号  sigx sigy  tau  sigx  sigy  tau   sigx   sigy    tau  sigx   sigy    tau \n')
    end
    Cxi=[ gx(1), gx(2), gx(2), gx(1)];                %工作数组，用于存计算位置点坐标
    Eta=[ gx(1) , gx(1), gx(2), gx(2)];
    Stress = zeros(ne,12);          %存储单元的4个高斯点或结点的主应力和Mises应力
    P_Stress = zeros(12,1);         %暂存单元的4个计算点的应力分量，便于输出
    ES=zeros(3,1);               % Es ---单元应力分量列阵（3×1）：[sigx,  sigy,  tauxy]
    De=zeros(8,1);               %单元结点位移向量
    D = Elastic_Matrix (pm, E, nv);    %调用函数的弹性矩阵，E-弹性模量，nv—泊松比
    for ie=1:1:ne
        for r=1:1:4                                %从结点向量中提取单元结点位移向量
            i0=2* EL(ie,r);
            De(2*r-1,1)=U(i0-1);
            De(2*r,1)=U(i0);
        end
        for r=1: 4
            [B, dt] =Quad_B4_Matrix( ie, Cxi (r), Eta (r) );     %调用函数计算几何矩阵、雅克比行列式
            S=D*B;
            ES=S*De ;                                    %单元的应力分量
            Sig12 = Main_Strees(ES(1,1),ES(2,1),ES(3,1));   %调用函数：计算单元主应力和Mises应力
            if  g_node ~=1
                j0=3*(r-1)+1;
                St_nd(ie,j0:(j0+2))=ES(1:3,1);
                %      St_ndP(ie,j0:(j0+2))= Sig12(1:3);
            end
            for s=1:1:3                              %将单元的应力分量转存在Stress前3列
                j0=3*(r-1)+s;
                P_Stress(j0,1)=ES(s,1);                  %将应力分量转存到显示用数组
                Stress(ie,j0)= Sig12(1,s);               %将主应力、Mises应力存到矩阵
            end
        end
        fprintf(fid,[repmat('%5i ', 1, 1) repmat('%12.4f ',1,12)],ie, P_Stress(:,1))    %12个应力分量
        fprintf(fid,' \n')
    end
    if g_node ==1
        fprintf(fid,'\n          单元高斯积分点主应力及Mises应力 \\n')
        fprintf(fid,' 1#（-0.57735，-0.57735） 2#（0.57735，-0.57735） 3#（0.57735，0.57735）4#（-0.57735，0.57735）\n')
        fprintf(fid,'单元号   Sig1   sig2  Mises  Sig1   sig2    Mises   Sig1  sig2   Mises  Sig1  sig2   Mises \n')
    else
        fprintf(fid,'\n           单元结点处主应力及Mises应力 \n')
        fprintf(fid,'\ n     1#（-1，-1）  2#（1，-1）     3#（1， 1）    4#（-1，1）\n');
        fprintf(fid,'结点号 Sig1   sig2  Mises  Sig1   sig2    Mises   Sig1  sig2   Mises  Sig1  sig2   Mises \n')
    end
    for ie=1:ne                            %输出所有主应力及Mises应力
        fprintf(fid,[repmat('%5i ', 1, 1) repmat('%12.4f ',1,12)],ie, Stress(ie,:))    %12个应力分量
        fprintf(fid,' \n')
    end
end

if  g_node ~=1
    %计算绕结点的应力，以下部分可调用函数Node_Stress(Stress，AA)实现
    fprintf(fid,'\n          绕结点应力计算结果 \n')  ;
    fprintf(fid,'Node    sigx        sigy         tau         sig1         sig3         Mises \n');
    Stress_nd = zeros(nd,6 ) ;
    for  i=1:1:nd
        Sd = zeros( 1, 3 ) ; num=0;
        for ie=1:1:ne
            for k=1:1:4
                if EL(ie,k) == i
                    j0=3*(k-1)+1;
                    Sd = Sd + St_nd(ie,j0:(j0+2));
                    num=num+1;
                end
            end
        end
        Stress_nd(i,1:3) = Sd / num ;
        Sig12 = Main_Strees(Stress_nd(i,1),Stress_nd(i,2),Stress_nd(i,3));
        Stress_nd(i,4:6) = Sig12(1:3);
        fprintf(fid,[repmat('%5i ', 1, 1) repmat('%12.4f ',1,6)],i,Stress_nd(i,:));
        fprintf(fid,' \n');
    end
    fprintf(fid,' \n  绕点平均应力结果统计分析 \n');
    s_Max=max(Stress_nd);
    fprintf(fid,[repmat('%s ', 1, 1) repmat('%12.4f ',1,6)],'最大值', s_Max);     %结点应力各个量的最大值
    fprintf(fid,' \n');
    s_Min=min(Stress_nd);
    fprintf(fid,[repmat('%s ', 1, 1) repmat('%12.4f ',1,6)],'最小值', s_Min);     %结点应力各个量的最小值
    fprintf(fid,' \n');
    
end
return


















