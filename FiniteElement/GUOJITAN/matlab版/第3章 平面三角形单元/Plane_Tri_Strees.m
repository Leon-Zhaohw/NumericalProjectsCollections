

 	function Stress_nd = Plane_Tri_Strees(file_out,U)
 	%  求单元应力分量及主应力、Mises应力
 	%  输入参数：U----总的结点位移向量； file_out—字符串，存储计算结果的文件名
 	% 返回量：绕点平均后的结点应力分量、主应力、Mises应力
 	 global pm E nv t0  lou nd ne ng nj    XY EL    BC  QJ    
 	%全局变量：结点数、单元数、结点坐标、单元信息 、弹性矩阵 
 	     fid=fopen(file_out,'at');                %打开并追加存储应力计算结果的文件   
 	    %计算单元应力
 	  Stress = zeros(ne,6);         
 	  AA=zeros(ne,1);              
 	   ES=zeros(3,1);              
 	   De=zeros(6,1);             
 	fprintf(fid,'\n          单元应力计算结果 \n') ;  
 	fprintf(fid,'\n Element    sigx        sigy         tau         sig1         sig3         Mises \n'); 
 	  D = Elastic_Matrix (pm, E, nv);
 	 for ie=1:1:ne
 	      for r=1:1:3                               
 	       i0=2* EL(ie,r);r2=2*r;
 	       De([r2-1,r2],1)=U([i0-1,i0],1); 
 	      end
 	
 	   [B, A] =  Plane_B3_Matrix( ie);           
 	   AA(ie,1)=A;
 	   S = D * B ;                                    
 	   ES=S*De  ;                      
 	   Stress(ie,1:3)=ES(1:3,1);                        
 	    Sig12 = Main_Strees(ES(1,1),ES(2,1),ES(3,1)); 
 	     Stress(ie,4:6) = Sig12(1:3);                      
 	   fprintf(fid,[repmat('%5d ', 1, 1) repmat('% 4f ',1,6)],ie, Stress(ie,:));  
 	  fprintf(fid,' \n');  
 	 end  
 	  fprintf(fid,' \n  应力结果统计分析 \n'); 
 	  s_Max=max(Stress);
 	  fprintf(fid,[repmat('%s ', 1, 1) repmat('% 4f ',1,6)],'最大值', s_Max);     %输出应力各个量的最大值
 	  fprintf(fid,' \n'); 
 	   s_Min=min(Stress);
 	  fprintf(fid,[repmat('%s ', 1, 1) repmat('% 4f ',1,6)],'最小值', s_Min);     %输出应力各个量的最小值
 	  fprintf(fid,' \n'); 
 	 
 	  %计算绕结点的应力，以下部分可调用函数Node_Stress(Stress，AA)实现
 	fprintf(fid,'\n          绕结点应力计算结果 \n')  ; 
 	fprintf(fid,'Node    sigx        sigy         tau         sig1         sig3         Mises \n'); 
 	 Stress_nd = zeros(nd,6 ) ;           
 	   for  i=1:1:nd                                
 	     Sd = zeros( 1, 3 ) ;                        
 	     A_tol = 0 ;
 	   for ie=1:1:ne
 	     for k=1:1:3
 	        if EL(ie,k) == i                    
 	           Sd = Sd + Stress(ie, 1:3 ) * AA(ie); 
 	           A_tol = A_tol + AA(ie);
 	        end
 	      end
 	    end
 	      Stress_nd(i,1:3) = Sd / A_tol ;             
 	      Sig12 = Main_Strees(Stress_nd(i,1),Stress_nd(i,2),Stress_nd(i,3));
 	      Stress_nd(i,4:6) = Sig12(1:3);                             
 	     fprintf(fid,[repmat('%5i ', 1, 1) repmat('% 4f ',1,6)],i,Stress_nd(i,:));
 	     fprintf(fid,' \n');
 	   end
 	      fprintf(fid,' \n  绕点平均应力结果统计分析 \n');
 	     s_Max=max(Stress_nd);
 	  fprintf(fid,[repmat('%s ', 1, 1) repmat('% 4f ',1,6)],'最大值', s_Max);     %结点应力各个量的最大值
 	  fprintf(fid,' \n'); 
 	  s_Min=min(Stress_nd);
 	  fprintf(fid,[repmat('%s ', 1, 1) repmat('% 4f ',1,6)],'最小值', s_Min);     %结点应力各个量的最小值
 	  fprintf(fid,' \n');
 	 
 	fclose all;
 	return









