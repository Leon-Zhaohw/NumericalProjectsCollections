%WRITE DATA TO FILE (for other purposes later)

function writeToFile(x1)

% Write data to file (as a string of floating point numbers) 
% Function is written by AIV. Version 15. October 2017

filename = input('Give the name of new file for storing ...
	results: ', 's');

fid = fopen(filename,'w');
fwrite(fid,x1(:,1),'double');
fwrite(fid,x1(:,2),'double');
status = fclose(fid);
return;