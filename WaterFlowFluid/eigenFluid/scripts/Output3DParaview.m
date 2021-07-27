[x,y,z] = meshgrid(0:pi/50:pi, 0:pi/50:pi, 0:pi/50:pi);
%Two Neumann wall on x axis and Two Neumann wall on y axis.

vx = cos(x).*sin(y).*cos(z);
vy = sin(x).*cos(y).*cos(z);
vz = 2*sin(x).*sin(y).*sin(z);

B = cat(4,vx,vy,vz);
B = permute(B, [4 1 2 3]);
fid = fopen('TwoNeumanXY.raw','w');
fwrite(fid,B,'float'); 
fclose(fid);