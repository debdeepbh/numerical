% Get optimal scaling values of wpx
function wpxsc = getwpxsc(type, p, sigma)

getdata;
wpxsc = zeros(8,p+1);

for i=1:8
	z = wpx(i,:);
	z = z(1:1024);

	wpxsc(i,:) = getoptsc(z, K, type, p, sigma, 'bisec');
	i
end


csvwrite(strcat('wpxsc_',type,num2str(p)), wpxsc);

figure;
for i=1:p+1
	subplot(p+1,1,i);
	hist(wpxsc(:,i));
	xlim([0 1]);
	title(num2str(i));
end
