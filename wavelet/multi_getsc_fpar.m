% get the optimal scaling values for multi_fone()

function alpha = multi_getsc_fpar(wax, aximp, testyori, type, p, noiseax)

errorl = 0.04;

rho = 1;
method = 'hard'; 	%for theoretical puroposes

alpha = zeros(1,p+1)+1; % start with scaling value 1

% use bisection method


for i=1:p+1
	% initial guess
	a = 0.1;
	b = 0.9;
	while (abs(a-b)>errorl)
	alpha(i) = a;
	[z, ratiounthr] = multi_fone(wax, aximp, testyori, type, p, noiseax, alpha, rho, method);
	fa = ratiounthr(i);

	alpha(i) = b;
	[z, ratiounthr] = multi_fone(wax, aximp, testyori, type, p, noiseax, alpha, rho, method);
	fb = ratiounthr(i);

	
	% if the signs are opposite
	if ((a-fa)>=0 && (b - fb)<=0)
		% disect
		c = (a+b)/2;
		alpha(i) = c;
	[z, ratiounthr] = multi_fone(wax, aximp, testyori, type, p, noiseax, alpha, rho, method);
		fc = ratiounthr(i);
		if ((c -fc) >0)
			a = c;
		else
			b = c;
		end
	% if the signs are the same
	elseif ((a-fa)<=0 && (b - fb)>=0)
		c = (a+b)/2;
		alpha(i) = c;
	[z, ratiounthr] = multi_fone(wax, aximp, testyori, type, p, noiseax, alpha, rho, method)
		fc = ratiounthr(i);
		if ((c -fc) >0)
			b = c;
		else
			a = c;
		end
	else
		i
		c = 0.1;
		disp('No solution in between');
		break;
	end
	end
	alpha(i) = c;
end


