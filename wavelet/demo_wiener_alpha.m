% 
j=1;
sigmalist = 1:1:50;
for sigma=sigmalist
	thenoise = randn([1 1024])*sigma;
	sig = testconv + thenoise;
A = 0.01:0.5:40;
Ns = length(A);
i=1;
errwi = zeros(1,Ns);
for alpha=A
	
	%subplot(Ns,1,i);
	dw = decwien(sig,K,thenoise,alpha);
	errwi(i) = rele(dw,testyori);
	%plot(dw);
	i=i+1;
end
plot(A, errwi)
hold on
minA(j) = A(find(errwi == min(errwi)))(1);
mine(j) = min(errwi);
j=j+1;
end

figure
plot(sigmalist,minA)
