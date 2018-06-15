function demo_save(fignum, filename, sizeval)
set(fignum,'PaperPositionMode','auto');  
pathname = '/home/debdeep/gdrive/anita/frankdata/pic/';
fullpath = strcat(pathname,filename,'.pdf');
print(fignum,'-dpdf',sizeval,fullpath)
