% AMG_USERSET_MAT sets us up for User Specified Matrix

amg_globals;

dir = input('\nPlease enter the directory name containing the binary crs matrix.\n:', 's');

fr = fopen(strcat(dir,'col.bin'),'r');
i = fread(fr,'int32');   %row
fclose(fr);

%fc = fopen(strcat(dir,'row.bin'),'rb');
%j = fread(fc,'int32');   %col
%fclose(fc);

fc = fopen(strcat(dir,'cnt.bin'),'r');
c = fread(fc,'int32');   %cnt
fclose(fc);

fe = fopen(strcat(dir,'ele.bin'),'r');
a = fread(fe,'double');   %ele
fclose(fe);

fb = fopen(strcat(dir,'b.bin'),'r');
b = fread(fb,'double');   %b
fclose(fb);

j = zeros(1,length(i));
p = 1; q = 1;
for k = 1 : length(c)
    for l = 1 : c(k)
        j(q) = p;
        q = q + 1;
    end
    p = p + 1;
end

i = i + 1;
%A = sparse(i,j,a);

FINEPOINTS = length(c);
A(1).matrix = sparse(i,j,a); 

X_Guess = zeros( length(c), 1 );

RHS = ones( length(c), 1 );
