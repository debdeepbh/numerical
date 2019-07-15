function  genparams(Nbd) 
% computes the parameters before running the time integration
% generate global variables

%% Load these from the file
load('Lam_00.mat')
load('Lam_10.mat')
load('Lam_20.mat')

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
