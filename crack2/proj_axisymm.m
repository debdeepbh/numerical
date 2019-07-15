function [out, u0_r, u0_z] = proj_axisymm(initial)

global length;
global width;
global dx;
global dy;
global delta;

%global Nbd;

Nbd = initial;

global Xdiff1;
global Xdiff2;
global Xdiffnorm;

global Lam_00;
global Lam_10;
global Lam_20;

global r
global r_p

global z
global z_p

%global z_diff
%global rrpzdiffsq	% no need
global neg2rrp 		% - 2 * r*r_p

nx = length/dx+1;
ny = width/dy+1;
totalnodes = nx*ny;
%ind =1;

%initial conditions

sigma=12e6/dx;
%sigma=24e6/dx;

rho=2440;
Gnot = 135;
E = 72e9;
% changed here
nu = 0.33;       % comment on paper: for 2d, the effective poisson ration is 1/3
%nu = 0.25;

snot = sqrt(4 * pi * Gnot /(9*E*delta));

%% cnot is the PMB constnat
cnot = 6*E/( pi * (delta^3) * (1 - nu));

%% Multiply the lambdas with the PMB constant
lam_00 = Lam_00 * cnot;
lam_10 = Lam_10 * cnot;
lam_20 = Lam_20 * cnot;

%f = figure('visible','off');

%internal force matrices
ss = size(Nbd);
m = ss(1);
n = ss(2);

serCol = (1:m)'; %  a list of all serials
Cent = repmat(serCol, 1, n);	% a matrix of all serials of the centers



%% External force goes here
extforce_r = zeros(m,1); % r component of  external force on i-th serial
extforce_z = zeros(m,1); 

%% External force is only in the z direction
%% Applies only on the single outer layer of the body

for j=1:nx
    extforce_z(ind2ser(1,j)) = -(sigma);
    extforce_z(ind2ser(ny,j)) = sigma;
end

uold_r = zeros(m,1); % this vector stores the r component of the old displacement of all points, to get the displacement of value of i-th element, use uold_r(i)
uold_z = zeros(m,1);

uolddot_r = zeros(m,1);
uolddot_z = zeros(m,1);
uolddotdot_r = zeros(m,1);
uolddotdot_z = zeros(m,1);




dt = 25e-9;
tic

%for t = 1:1800
for t= 1:2000


    % loop starts
    u0_r = uold_r + dt * uolddot_r + dt * dt * 0.5 * uolddotdot_r;
    u0_z = uold_z + dt * uolddot_z + dt * dt * 0.5 * uolddotdot_z;

    % u0 is the long vector with `totalnodes` rows and 2 columns, stores the displacements of all nodes
    
    mask = (~~Nbd);


     
    %for i = 1:totalnodes	% i is the serial of the center

    %    % UI is the displacement of the center
    %    u_r = u0(i,1);
    %    u_z = u0(i,2);

    % dynamically computed displacement vectors of the centers
    u_r = u0_r(Cent);
    u_z = u0_z(Cent);

    % dynamically computed displacement vectors of the neighbors 
    u_r_p = u0_r(Nbd + (~mask));
    u_z_p = u0_z(Nbd + (~mask));


    %% Define a few frequently used quantities
    zdiff = z_p - z;

    uzdiff = u_z_p  - u_z;

    xietaz = zdiff .* uzdiff;	% xi_z eta_z

    ruord = r_p.*u_r_p + r.*u_r;	% r u_r + r' u_{r'}
    rurev = r .* u_r_p + r_p .* u_r;	% r'u_r + r u_{r'}

    bigrterms = (r.^2 + r_p.^2) .* u_r_p  + (-neg2rrp) .* u_r;



    %% The integrands without the r' from the Jacobian

    intforce_r = (-r).* ( ruord + xietaz) .* lam_00  + (bigrterms + r_p .* xietaz) .* lam_10 - r_p .* rurev .* lam_20;

    intforce_z = zdiff .* ( ruord + xietaz) .* lam_00 + zdiff .* rurev .* lam_10;

    % Multiply by the Jacobian part: r'
    intforce_r = intforce_r .* r_p;
    intforce_z = intforce_z .* r_p;

    % zero out the missing values
    intforce_r = intforce_r .* mask; 
    intforce_z = intforce_z .* mask; 

    % add all the force contributions from all neighbors
    % adding all the values in the same row
    totalintforce_r = sum(intforce_r,2);
    totalintforce_z = sum(intforce_z,2);

    %    totalintforce(i,1) = 0; % this is where I will store my total internal force
    %    totalintforce(i,2) = 0; % this is where I will store my total internal force
    %    
    %    %% debug
    %    i

    %    for j=1:n		
    %        serJ = Nbd(i,j);	% serJ is the serial of a generic neighbor
    %        if serJ > 0	% only existing bonds
    %    	
    %    	    %%  Position of the point j
    %    	    posJ = ser2pos(serJ);
    %    	    r_p = posJ(1);
    %    	    z_p = posJ(2);

    %    	    % UJ is the displacement of the point j
    %    	    u_r_p = u0(serJ,1);
    %    	    u_z_p = u0(serJ,2);

    %    	    %% lambdas
    %    	    %[lam_00, lam_10, lam_20] = params(i, serJ);
    %    	    lam_00 = Lam_00(i,j);
    %    	    lam_10 = Lam_10(i,j);
    %    	    lam_20 = Lam_20(i,j);

    %    	    intforce_r = r_p * ( ((-1) *r) * (r_p * u_r_p + r* u_r + (u_z_p - u_z) * (z_p - z) ) * lam_00 + (u_r_p *(r^2 + r_p^2) + 2 * r * r_p * u_r + r_p * (z_p - z) * (u_z_p - u_z)) * lam_10 - (r * r_p * u_r_p + r_p^2 * u_r) * lam_20 );
    %    	    % need to add body force now

    %    	    intforce_z = r_p *( (z_p - z) * ( r_p* u_r_p + r* u_r + (z_p - z)* (u_z_p - u_z)) * lam_00 - (z_p -z)* (r * u_r_p + r_p * u_r) * lam_10 );
    %    	    % need to add body force now

    %    	    
    %    	    totalintforce(i,1) = totalintforce(i,1) +intforce_r;
    %    	    totalintforce(i,2) = totalintforce(i,2) +intforce_z;


    %    	    %%%%%%% debug
    %           % Udiff1(i,j) = u0(Nbd(i,j),1) - u0(i,1);
    %           % Udiff2(i,j) = u0(Nbd(i,j),2) - u0(i,2);
    %        end
    %    end
    %end
    
%    Bigvect1 = Udiff1 + Xdiff1;
%    Bigvect2 = Udiff2 + Xdiff2;
%    
%    Bigvectnorm = sqrt(Bigvect1.*Bigvect1 + Bigvect2.*Bigvect2);    % norm(ucap - u + xacp -x)
%    S = (Bigvectnorm - Xdiffnorm ) ./ Xdiffnorm;    %gives NaN when an element  is not a nbd
%    
%    Multipl = cnot * (S ./ Bigvectnorm);
%    Force1 = Multipl .* Bigvect1;
%    Force2 = Multipl .* Bigvect2;
    
%    totalintforce = [sum(Force1.*restrict, 2) sum(Force2.*restrict,2)];    %sum(A) ad totalintforce2 = sum(Force2, 2);
    
    u0dotdot_r = (1/rho) .*( totalintforce_r .* (dx.*dy) + extforce_r) ;
    u0dotdot_z = (1/rho) .*( totalintforce_z .* (dx.*dy) + extforce_z) ;

    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Bond breaking criteria
    %%%%%%%% disabled %%%%%%
    % Nbd = ((S - snot ) < 0).* Nbd;
    %%%%%%%%%%%%%%%%%%%%%%5
    
    %calculate u0dot
    u0dot_r = uolddot_r + dt * 0.5 * uolddotdot_r + dt * 0.5 * u0dotdot_r;
    u0dot_z = uolddot_z + dt * 0.5 * uolddotdot_z + dt * 0.5 * u0dotdot_z;
    
    uold_r = u0_r;
    uold_z = u0_z;
    uolddot_r = u0dot_r;
    uolddot_z = u0dot_z;
    uolddotdot_r = u0dotdot_r;
    uolddotdot_z = u0dotdot_z;
    
    %loop ends
    t
%     if mod(t,50) == 0
%         drawmesh(Nbd);
%         M(ind) = getframe(f);
%         ind = ind +1;
%     end
end
toc
out = Nbd;
%drawmesh(out);
%movarr = M;



