[x,y,z] = meshgrid(0:pi/50:pi, 0:pi/50:pi, 0:pi/50:pi);
%Two Neumann wall on x axis and Two Neumann wall on y axis.

vx = sin(x).*cos(y).*cos(z);
vy = cos(x).*sin(y).*cos(z);
vz = -2*cos(x).*cos(y).*sin(z);

B = cat(4,vx,vy,vz);
B = permute(B, [4 1 2 3]);
fid = fopen('Dirichlet.raw','w');
fwrite(fid,B,'float'); 
fclose(fid);