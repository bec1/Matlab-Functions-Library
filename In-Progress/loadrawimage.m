function num=loadrawimage(filename)
%read the absorption image from aia or fits file
% For aia files, load the three base images (with atoms, without atoms,
% background), calculate the absorption image and prepare everything
% for display
% For fits files, just load the image

sourcepath = ('\\Elder-pc\j\Elder Backup Raw Images');
year = filename(7:10);
month = filename(1:2);
day = filename(4:5);
if str2num(filename(12:13))<7 %if we're doing super late night stuff, we've slipped back in time
    day = num2str(str2num(day)-1,'%02.f');
end
yrmnth = strcat(year,'-',month);
yrmnthday = strcat(yrmnth,'-',day);
filepath = strcat(sourcepath,'\',year,'\',yrmnth,'\',yrmnthday);
filename = strcat(filepath,'\',filename,'.fits');


img=fitsread(filename);
pxsize = 1;
num = AtomNumber(img,pxsize^2,0.215,inf);
     
end