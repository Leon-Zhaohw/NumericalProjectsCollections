

 
	function [file_in , file_out] = File_Name
	  %数据文件管理，输入原始数据文件路径及文件名，并在同路径下自动生成或自行定义输出文件名
	  %检查到自动生成的输出文件名已存在，则提醒覆盖或新建文件名，也可修改文件路径。
	  [filename, path_str] = uigetfile( '../*.xls;*.xlsx;*.txt',' 选择有限元模型数据文件')
    [s,name_str,ext_str] = fileparts( filename );
	    file_in =fullfile( path_str, filename ) ;
	    ext_str_out = '.txt' ; 
	  file_out = fullfile( path_str, [name_str, ext_str_out] )     
	      % 检查输出文件是否存在
   while  exist( file_out ) ~= 0 
	      ss=strcat('拟保存计算结果的文件名：',file_out,' 已存在，是否覆盖？');
       button=questdlg(ss,'是否覆盖已有文件？','覆盖','新建','新建');
     if button == '新建'
	          [name_str,path_str]=uiputfile( '../*.txt',' 保存有限元计算结果文件名')
        file_out = fullfile( path_str,name_str )
     elseif button == '覆盖'
        return
    end 
 end
return
