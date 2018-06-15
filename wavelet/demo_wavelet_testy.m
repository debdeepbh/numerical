plotcoeffs(wtrans(testy(1:1024),'d10',4),4);

h = gcf;
set(h,'PaperPositionMode','auto');  
print gcf -dpdf '-S800,1200' '/home/debdeep/gdrive/anita/frankdata/pic/demo_wavelet_testy.pdf'
