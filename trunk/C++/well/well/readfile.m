
fidp = fopen('Efield2_299.jd','r','l');
if fidp==-1
    a=1;
    return
end

var=fread(fidp,1000,'double')
plot(var);
fclose(fidp);