function str = txtread(filename)

fid = fopen(filename); 
str = fread(fid,[1,inf],'char'); 
fclose(fid);
str = char(str);

