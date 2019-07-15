function rndtest(Nbd)

ss = size(Nbd);

dummy = zeros(ss);


center1 = ser2pos1(Nbd);
generic1= ser2pos1((1:ss(1))');
absval = abs(center1 - generic1);

absval = absval.*(~~Nbd);
max(max(absval))

maxnorm = 0;
maxabs = 0;
for r=1:ss(1)
%r = floor(rand()*10000)
for i=1:ss(2) 
	s = Nbd(r,i);
	if s > 0
		dummy(r,i) =1;
		%maxnorm = max(norm(ser2pos(r) - ser2pos(s)), maxnorm);
		currnorm = sqrt((ser2pos1(r) - ser2pos1(s))^2 + (ser2pos2(r) - ser2pos2(s))^2);
		currabs =  abs(ser2pos1(r) - ser2pos1(s));

		maxnorm = max(currnorm, maxnorm);
		maxabs = max(currabs, maxabs);
		if maxabs > 0.002
			currnorm
			r
			i
			s
			return
		end
	end
end
end

maxabs

sum(sum((dummy ~= ~~Nbd)))
ss(1)*ss(2)
