% plot the histogram of optimal scaling parameters
for i=1:p+1
	subplot(p+1,1,i)
	hist(sc_wpx(:,i));
end

