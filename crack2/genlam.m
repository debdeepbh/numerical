function  genlam(Nbd) 
% computes the parameters between the points with serial i and j



% this is where we store the Lambdas
global Lam_00
global Lam_10
global Lam_20

global r
global z

global r_p
global z_p

global z_diff
%global rrpzdiffsq	% no need
global neg2rrp



global delta


mask = (~~Nbd);	% this matrix has zero in the missing positions

ss = size(Nbd);	% get the size of the nbd array

Lam_00 = zeros(ss);
Lam_10 = zeros(ss);
Lam_20 = zeros(ss);

serList = (1:ss(1))';

r = repmat(ser2pos1(serList), 1, ss(2)).*mask;
z = repmat(ser2pos2(serList), 1, ss(2)).*mask;


%r = repmat(PosVec1,1,ss(2)); % repeat the column vector and expand it to the size of Nbd
%z = repmat(PosVec2,1,ss(2)); % repeat the column vector and expand it to the size of Nbd

r_p =  ser2pos1(Nbd).*mask;
z_p =  ser2pos2(Nbd).*mask;


z_diff = z_p - z;

%%% Checked that the differences are within delta
max(max(z_diff))
min(min(z_diff))

r_diff = abs(r_p - r);

min(min(r_diff))
max(max(r_diff))
%%[row, col] = find(r_diff > delta)



rrpzdiffsq = r.^2 + r_p.^2+ z_diff.^2;
neg2rrp = (-2).* r .* r_p;

% to compute K, we are dividing by b, so replace missing values by 1
neg2rrp = neg2rrp.*mask + (~mask);
% Note that, b can still be zero if one of r or r_p is zero. This is why 
% we need to consider a problem where all the points are away from the axis of symmetry

%% Checking that b does not produce zero values
min(min(abs(neg2rrp)))

%for i = 1:ss(1)
%	for j =1:ss(2)
%		if Nbd(i,j) > 0
%			if b(i,j) == 0 
%				b(i,j)
%				i,j
%				r(i,j)
%				r_p(i,j)
%				return
%			end
%		end
%	end
%end






%K = (delta^2 - z_diff.^2 - r_p.^2 - r.^2)./((-2)*r.*r_p);
K = (delta^2 - rrpzdiffsq)./(neg2rrp);
K = K.* mask;

%% Checking that the values are between -1 and 1, so that we can compute
%% arccos(K) later
%% In particular, when away from the axis of symmetry, between 0 and 1
max(max(K))
min(min(K))

%% We shall replace the missing values by 1 so that acos will produce zero
K = K.*mask + (~mask);

%% Now compute the arccos
alpha = acos(K);

%% Check that the values  are between 0 and pi/2, when away from the axis of symmetry
min(min(alpha))
max(max(alpha))


% influence function related parameters 
fun_00 = @(phi, a, b) ( 1./ (( a +  b .* cos(phi) ).^(3/2)) );
fun_01 = @(phi, a, b) ( cos(phi)./ ((a + b .* cos(phi) ).^(3/2)) );
fun_20 = @(phi, a ,b) ( (cos(phi).^2)./ ( ( a + b .* cos(phi)).^(3/2)) );


% Does the integrand and limit both have the same expression inside?

% Do we _have_ to do a term-wise integration?
tic
for i = 1:ss(1)

	% print i
	if ~mod(i,10)
		i
	end

	for j = 1:ss(2)
		if Nbd(i,j) > 0
			% spot the bad apples
			%if K(i,j) >= 1
			%	j
			%end
	%Lam_00 = 2 *integral(@(phi)fun_00(phi, rrpzdiffsq(i,j), neg2rrp(i,j)), 0, alpha(i,j));
	%Lam_10 = 2 *integral(@(phi)fun_01(phi, rrpzdiffsq(i,j), neg2rrp(i,j)), 0, alpha(i,j));
	%Lam_20 = 2 *integral(@(phi)fun_20(phi, rrpzdiffsq(i,j), neg2rrp(i,j)), 0, alpha(i,j));

	Lam_00(i,j) = integral(@(phi)fun_00(phi, rrpzdiffsq(i,j), neg2rrp(i,j)), (-1)*alpha(i,j), alpha(i,j));
	Lam_10(i,j) = integral(@(phi)fun_01(phi, rrpzdiffsq(i,j), neg2rrp(i,j)), (-1)*alpha(i,j), alpha(i,j));
	Lam_20(i,j) = integral(@(phi)fun_20(phi, rrpzdiffsq(i,j), neg2rrp(i,j)), (-1)*alpha(i,j), alpha(i,j));
		end
	end
end
toc

save('Lam_00.mat','Lam_00')
save('Lam_10.mat','Lam_10')
save('Lam_20.mat','Lam_20')

%Lam_00 = integral(@(phi)fun_00(phi, a, b), (-1)*alpha, alpha);
%Lam_10 = integral(@(phi)fun_01(phi, a, b), (-1)*alpha, alpha);
%Lam_20 = integral(@(phi)fun_20(phi, a, b), (-1)*alpha, alpha);
