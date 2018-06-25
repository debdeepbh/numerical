% get optimal scaling for different noise variance
csig = testconv;

smax = 50;
p = 4;

x = [testvec zeros(1,512)];
nx = norm(x);

sc = zeros(smax,p+1);

for sigma=1:smax
	ns = randn([1 1024])*sigma; 
	sig = csig + ns;
	sc(sigma,:) = getoptsc(sig, K, type, p, ns, 'bisec')
end

for i=1:p+1
	subplot(p+1,1,i);
	hist(sc(:,i));
	title(num2str(i));
end

