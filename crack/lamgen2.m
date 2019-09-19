function lamgen2(Nbd,dx, dy, nx, ny, delta)


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

Lam_00 = zeros([m,n]);
Lam_10 = zeros([m,n]);
Lam_20 = zeros([m,n]);


%% Define functions, algebraic
fun_00 = @(beta, pos2rrp, extra)( 4./ ((beta.^2)  .* sqrt(pos2rrp.^2  - (extra - beta.^2))));
fun_10 = @(beta, pos2rrp, extra)( 4 *( extra - beta.^2)./ ((pos2rrp).*(beta.^2) .* sqrt(pos2rrp.^2  - (extra - beta.^2))));
fun_20 = @(beta, pos2rrp, extra)( (4 *( extra - beta.^2).^2) ./ ( ((pos2rrp).^2).*(beta.^2)  .* sqrt(pos2rrp.^2  - (extra - beta.^2))));

%% Do the integration, term-wise
tic
for i=1:m
	i
	for j = 1:n
		if Nbd(i,j) > 0
			%if j ==1
			%	xi_norm(i,j)
			%	delta
			%	return
			%end
			Lam_00(i,j) = integral(@(beta)fun_00(beta, pos2rrp(i,j), extra(i,j)), xi_norm(i,j), delta);
			Lam_10(i,j) = integral(@(beta)fun_10(beta, pos2rrp(i,j), extra(i,j)), xi_norm(i,j), delta);
			Lam_20(i,j) = integral(@(beta)fun_20(beta, pos2rrp(i,j), extra(i,j)), xi_norm(i,j), delta);
		end
	end
end
toc

save('Lam2_00.mat','Lam_00')
save('Lam2_10.mat','Lam_10')
save('Lam2_20.mat','Lam_20')
%
%
