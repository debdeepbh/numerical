function lamgen(Nbd,dx, dy, nx, ny, delta)


[m,n] = size(Nbd);

%%%%%%%%%%%%%% Pre computation %%%%%%%

longvec = (1:m)';
r = repmat(ser2pos1(longvec, dx, nx), 1, n);
z = repmat(ser2pos2(longvec, dy, nx), 1, n);

r_p = ser2pos1(Nbd, dx, nx);
z_p = ser2pos2(Nbd, dy, nx);


xi_r = r_p - r;
xi_z = z_p - z;

%max(max(abs(xi_r .* ~~Nbd)))
%max(max(abs(xi_z .* ~~Nbd)))
%max(max(abs(sqrt(xi_r.^2 + xi_z.^2) .* ~~Nbd)))
%delta

pos2rrp = 2* r.* r_p;

rrpxi_sq = r.^2 + r_p.^2 + xi_z.^2;


xi_norm = sqrt(xi_r.^2 + xi_z.^2);
extra = r.^2  + r_p.^2 + xi_z.^2;

%% alpha matrix
% Make sure we are not dividing by zero
K = (rrpxi_sq - delta^2)./(pos2rrp + (~Nbd));


% mask K
K = K .* (~~Nbd);


%% Issue while using closed ball B_\delta(x)
%The only values in K.*(~~K) that will have value greater than 1 
% are the ones on the boundary, due to precision error. Probably 
% we can reset them to 1
% Check: values bigger than 1 indeed occur on the boundary of the ball
disp('before max, min')
max(max(K))
min(min(K))

% find where K > 1, after masking
%issue = K > 1;
%% when the values are close to 1 
issue = abs(K - 1) < 1e-10;
% zero out the problematic places, then fill them up with 1
K = K .* (~issue) + ~~issue;

disp('after max, min')
max(max(K))
min(min(K))
%return

alpha = acos(K); 

Lam_00 = zeros([m,n]);
Lam_10 = zeros([m,n]);
Lam_20 = zeros([m,n]);



%%% Define functions, trigonometric
fun_00 = @(phi, rrpxi_sq, pos2rrp)(2./( (rrpxi_sq - pos2rrp.* cos(phi)).^(3/2)) ); 
fun_10 = @(phi, rrpxi_sq, pos2rrp)((2 .* cos(phi))./( (rrpxi_sq - pos2rrp.* cos(phi)).^(3/2)) ); 
fun_20 = @(phi, rrpxi_sq, pos2rrp)((2 .* (cos(phi).^2))./( (rrpxi_sq - pos2rrp.* cos(phi)).^(3/2)) ); 

%%% Do the integration, term-wise
tic
for i=1:m
	fprintf(repmat('\b',1,30));
	fprintf('i = %d, perc = %d %%',i, floor(i/m*100));
	for j = 1:n
		if Nbd(i,j) > 0

			Lam_00(i,j) = integral(@(phi)fun_00(phi, rrpxi_sq(i,j), pos2rrp(i,j)), 0, alpha(i,j));
			%disp('in')
			%if j ==1
			%	Nbd(i,j)
			%	alpha(i,j)
			%	Lam_00(i,j)
			%	return
			%end
			Lam_10(i,j) = integral(@(phi)fun_10(phi, rrpxi_sq(i,j), pos2rrp(i,j)), 0, alpha(i,j));
			Lam_20(i,j) = integral(@(phi)fun_20(phi, rrpxi_sq(i,j), pos2rrp(i,j)), 0, alpha(i,j));
		end
	end
end

toc
save('Lam_00_prenotch.mat','Lam_00')
save('Lam_10_prenotch.mat','Lam_10')
save('Lam_20_prenotch.mat','Lam_20')



