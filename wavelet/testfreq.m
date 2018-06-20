type = 'shan';
level = 2;
M = getbasismat(type, level,1024);

p = proj(z,type, level);

plotfanita(z-p);
% hold on here
hold on

plotfanita(p);
legend('z-p','p')


for j=1:level+1
	plotfanita(1000*M(j,:));
end


title(strcat('projection p of signal z onto *', num2str(level), '* stage *', type, '* wavelet space'))
hold off
