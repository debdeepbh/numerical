clear all
snrval_par = zeros(1,23);
snrval_one = zeros(1,23);
snrval_avg = zeros(1,23);
for indexthis=1:23
	[ax, aximp] = prepsig(indexthis);
	wienforwdant
	snrval_par(indexthis) = getsnr(z)
	snrval_one(indexthis) = getsnr(zzone)
	snrval_avg(indexthis) = getsnr(z_avg)
end


