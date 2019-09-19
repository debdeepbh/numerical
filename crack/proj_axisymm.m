function [out, u0] = proj_axisymm(initial, dx, dy, nx, ny, delta)

%%%%%%%%%%%%%%%%%%%%%
%%% Load these from the file
%% For Lam
load('Lam_00.mat')
load('Lam_10.mat')
load('Lam_20.mat')
%
%% For Lam2
%load('Lam2_00.mat')
%load('Lam2_10.mat')
%load('Lam2_20.mat')

% For Pre-notch
%load('Lam_00_prenotch.mat')
%load('Lam_10_prenotch.mat')
%load('Lam_20_prenotch.mat')

%%%%%%%%%%%%%%%%%%%%%%

imgcounter = 1;

Nbd = initial;

totalnodes = nx*ny;


%initial condition
%sigma=12e6/dx;
%sigma=12e6/dx;
%rho=2440;
%Gnot = 135;
%E = 72e9;
%nu = 0.22;       % commnet on paper: the effective poisson ration is 1/3
%
%snot = sqrt(4 * pi * Gnot /(9*E*delta));
%cnot = 6*E/( pi * (delta^3) * (1 - nu));

%% askari-silling paper, cylinder
nu = 0.25;
%rho = 2200e3;
rho = 2200;
K = 14.9e9;

%sigma=2e7/(2*pi*(0.1 + dx/2) * dx * dx );
sigma=2e6;


% From constant matching of PMB and classical strain energy
cnot = 18*K/pi/(delta^4);



%%multiply lambda by cnot, the PMB constant
lam_00= Lam_00 .* cnot;
lam_10= Lam_10 .* cnot;
lam_20= Lam_20 .* cnot;

%f = figure('visible','off');
f = figure('visible', 'off');

[m, n] = size(Nbd);

%%% Prescribed initial condition
maxdisp = 1e-4;
spread = 5;
%% A cos function with max height h at 1 and zero at L
cosdecay = @(x,h,L)( h/2.* cos(2*pi*(x-1)/2/(L-1)) + h/2);


inidisp = zeros(totalnodes,2); % extforce(i,:) = external force on i-th serial
extforce = zeros(totalnodes,2); % extforce(i,:) = external force on i-th serial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Force  %%%%%%%%%%%%%%%%%%%%%%

%%% Since our y-axis extends downward, the external forces directions are flipped in the $y$ direction.

%%% vertical force
%for j=1:nx
%    extforce(ind2ser(1,j,nx, ny),2) = -(sigma);	% force directions have to flipped
%
%    extforce(ind2ser(2,j,nx, ny),2) = -(sigma);	
%    extforce(ind2ser(ny,j,nx, ny),2) = sigma;
%
%    extforce(ind2ser(ny-1,j,nx, ny),2) = sigma;
%end

%%%% horizontal force, only on one layer
%for i=1:ny
%    extforce(ind2ser(i,1,nx, ny),1) = +(sigma);
%
%    %extforce(ind2ser(2,j,nx, ny),2) = -(sigma);	
%    %extforce(ind2ser(ny,j,nx, ny),2) = sigma;
%
%    %extforce(ind2ser(ny-1,j,nx, ny),2) = sigma;
%end

%%% horizontal force on several layers, reversed
%for i=1:ny
%	for lev = 1:spread
%	    dispval = cosdecay(lev, sigma, spread);
%
%	    %% reverse order in space to have zero force on the boundary
%	    extforce(ind2ser(i, spread - lev + 1 ,nx, ny),1) =  dispval;
%	end
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Displacements  %%%%%%%%%%%%%%

%%%% vertical displacement
%for j=1:nx
%	for lev = 1:spread
%	    dispval = cosdecay(lev, maxdisp, spread);
%
%	    inidisp(ind2ser(1 + lev -1,j,nx, ny),2) =  -dispval; % top layers, direction flipped
%	    inidisp(ind2ser(ny -lev +1,j,nx, ny),2) = dispval;
%    end
%end

%%% horizontal displacement
%for i=1:ny
%	for lev = 1:spread
%	    dispval = cosdecay(lev, maxdisp, spread);
%
%	    inidisp(ind2ser(i, 1 + lev -1,nx, ny),1) =  dispval;
%	    inidisp(ind2ser(i, nx -lev +1,nx, ny),1) = -dispval;
%	end
%end
%%% horizontal displacement, zero on the boundary
for i=1:ny
	for lev = 1:spread
	    dispval = cosdecay(lev, maxdisp, spread);

	    inidisp(ind2ser(i, spread - lev + 1,nx, ny),1) =  dispval;
	    %inidisp(ind2ser(i, nx -lev +1,nx, ny),1) = -dispval;
	end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





uold = zeros(totalnodes,2);
uolddot = zeros(totalnodes,2);
uolddotdot = zeros(totalnodes,2);

%%%%%%%%%%%%%% Pre computation %%%%%%%

longvec = (1:m)';
r = repmat(ser2pos1(longvec, dx, nx), 1, n);
z = repmat(ser2pos2(longvec, dy, nx), 1, n);

r_p = ser2pos1(Nbd, dx, nx);
z_p = ser2pos2(Nbd, dy, nx);

xi_r = r_p - r;
xi_z = z_p - z;

xi_norm = sqrt(xi_r.^2 + xi_z.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	Classical 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lame constants
%% Found from constant matching of linearized PMB and classical in cartesian coordinates
% Alternatively, mu = 3*K *(nu/(1+nu));
% which is same as below when nu=1/4;
mu = pi * delta^4 * cnot/30;

lamb = mu;

arow = zeros(1,nx);
acol = zeros(ny,1);

iMat = repmat((1:ny)', 1, nx);
jMat = repmat((1:nx), ny, 1);

serMat =  nx*(iMat - 1) + jMat;	% A matrix whose (i,j)th element is ind2ser(i,j)

% Use like displacement_r(serMat) to get u_r in a matrix format
u0Class_r = zeros(ny, nx);	% u_r values, in a matrix
u0Class_z = zeros(ny, nx);	% u_z values, in a matrix

rClass = r(serMat);
zClass = z(serMat);

% Required for time integration
uoldClass_r = zeros(ny, nx);
uoldClass_z = zeros(ny, nx);
uolddotClass_r = zeros(ny, nx);
uolddotClass_z = zeros(ny, nx);
uolddotdotClass_r = zeros(ny, nx);
uolddotdotClass_z = zeros(ny, nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 25e-9;
tic
for t = 1:3800

    % loop starts
    u0 = uold + dt * uolddot + dt * dt * 0.5 * uolddotdot;

    %%% Initial prescribed displacement
    if t==1
            u0 = inidisp;
    end

    %%% External force stops after a few time steps
    %if t==5
    %    extforce = zeros(totalnodes,2); % extforce(i,:) = external force on i-th serial
    %end

    %%% External force decreases in time, following a cosine curve
    %damptime =  100;
    %if t > damptime
    %        extforce = zeros(totalnodes,2);
    %else
    %        force_scaling = cosdecay(t, 1, damptime); 
    %        extforce = extforce .* force_scaling;
    %end

    
    restrict = (Nbd > 0);
    
    %%%%%%% New version
    disp_r = u0(:,1);
    disp_z = u0(:,2);

    %% Displacements of the centers
    u_r = disp_r(longvec);
    u_z = disp_z(longvec);

    %% Displacements of the neighbors
    u_r_p = disp_r(Nbd + (~Nbd));
    u_z_p = disp_z(Nbd + (~Nbd));

    eta_r = u_r_p - u_r;
    eta_z = u_z_p - u_z;


    uord = u_r .* r + u_r_p .* r_p;
    urev = u_r .* r_p + u_r_p .* r;

    uord_p_third = uord + eta_z .* xi_z;

    %% integrands of u_r and u_z without the jacobian
    int_r = uord_p_third .* (  r_p .* lam_10 - r .* lam_00 ) + urev .* ( r .* lam_10 - r_p .* lam_20);

    int_z = xi_z .* (uord_p_third .* lam_00 - urev .* lam_10);

    %% multiply by the Jacobian r_p
    Force1 = int_r .* r_p;
    Force2 = int_z .* r_p;

    %% mask it
    Force1 = Force1 .* (~~Nbd);
    Force2 = Force2 .* (~~Nbd);

    %%%% After masking
    
    totalintforce = [sum(Force1.*restrict, 2) sum(Force2.*restrict,2)];    %sum(A) ad totalintforce2 = sum(Force2, 2);
    
    u0dotdot = (1/rho) *( totalintforce * (dx*dy) + extforce) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 	Classical
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Borrowing prescribed external force
    extfr = extforce(:,1);
    extfz = extforce(:,2);

    extforce_r = extfr(serMat);
    extforce_z = extfz(serMat);


    %% These are matrices with dimension nx X ny
    u0Class_r = uoldClass_r + dt * uolddotClass_r + dt * dt * 0.5 * uolddotdotClass_r;
    u0Class_z = uoldClass_z + dt * uolddotClass_z + dt * dt * 0.5 * uolddotdotClass_z;


    %%% Add initial displacement
    if t==1
	    u0Class_r = reshape(inidisp(:,1), nx, ny)';
	    u0Class_z = reshape(inidisp(:,2), nx, ny)';
    end


    %% The last subscrips are for r and z, respectively
    % example_
    % ur_p0 is u_r(r_{j+1}, z_i)
    % ur_0n is u_r(r_{j}, z_{i-1})

    %% for u_r
    ur_00 = u0Class_r; 	% u_r(r_j, z_i)
    ur_p0 = shiftMat(u0Class_r, 0, 1, nx, ny);	% u_r(r_{j+1}, z_{i})		
    ur_n0 = shiftMat(u0Class_r, 0, -1, nx, ny);	% u_r(r_{j-1}, z_{i})		
    ur_n0 = shiftMat(u0Class_r, 0, -1, nx, ny);	% u_r(r_{j-1}, z_{i})		

    ur_0p = shiftMat(u0Class_r, 1, 0, nx, ny);	% u_r(r_{j}, z_{i+1})		
    ur_0n = shiftMat(u0Class_r, -1, 0, nx, ny);	% u_r(r_{j}, z_{i-1})		

    ur_pp = shiftMat(u0Class_r, 1, 1, nx, ny);	% u_r(r_{j}, z_{i+1})		
    ur_pn = shiftMat(u0Class_r, -1, 1, nx, ny);	% u_r(r_{j-1}, z_{i})		
    ur_np = shiftMat(u0Class_r, 1, -1, nx, ny);	% u_r(r_{j-1}, z_{i})		
    ur_nn = shiftMat(u0Class_r, -1, -1, nx, ny);	% u_r(r_{j-1}, z_{i})		

    %% for u_z
    uz_00 = u0Class_z; 	% u_r(r_j, z_i)
    uz_p0 = shiftMat(u0Class_z, 0, 1, nx, ny);	% u_r(r_{j+1}, z_{i})		
    uz_n0 = shiftMat(u0Class_z, 0, -1, nx, ny);	% u_r(r_{j-1}, z_{i})		

    uz_0p = shiftMat(u0Class_z, 1, 0, nx, ny);	% u_r(r_{j}, z_{i+1})		
    uz_0n = shiftMat(u0Class_z, -1, 0, nx, ny);	% u_r(r_{j}, z_{i-1})		

    uz_pp = shiftMat(u0Class_z, 1, 1, nx, ny);	% u_r(r_{j}, z_{i+1})		
    uz_pn = shiftMat(u0Class_z, -1, 1, nx, ny);	% u_r(r_{j-1}, z_{i})		
    uz_np = shiftMat(u0Class_z, 1, -1, nx, ny);	% u_r(r_{j-1}, z_{i})		
    uz_nn = shiftMat(u0Class_z, -1, -1, nx, ny);	% u_r(r_{j-1}, z_{i})		

    %% All the derivatives, finite difference
    %% dx is dr, dy is dz
    diff_r_ur = (ur_p0 - ur_00)/dx;
    diff_z_ur = (ur_0p - ur_00)/dy;

    diff_r_uz = (uz_p0 - uz_00)/dx;
    diff_z_uz = (uz_0p - uz_00)/dy;

    diff_rr_ur = (ur_p0 + ur_n0 - 2*ur_00) /(dx)^2;
    diff_zz_ur = (ur_0p + ur_0n - 2*ur_00) /(dy)^2;

    diff_rr_uz = (uz_p0 + uz_n0 - 2*uz_00) /(dx)^2;
    diff_zz_uz = (uz_0p + uz_0n - 2*uz_00) /(dy)^2;

    diff_rz_ur = (ur_pp - ur_pn - ur_np + ur_nn)/ (4 * dx *dy);
    diff_rz_uz = (uz_pp - uz_pn - uz_np + uz_nn)/ (4 * dx *dy);

    %% Forces
    forceClass_r = (lamb + mu) * ( (rClass .* diff_rr_ur + 2 * diff_r_ur)./ rClass - (rClass .* diff_r_ur + ur_00)./(rClass.^2) + diff_rz_uz) + mu * ( (rClass .* diff_rr_ur + diff_r_ur)./rClass - ur_00./(rClass.^2) + diff_zz_ur);

    forceClass_z = (lamb + mu) * ( (rClass .* diff_rz_ur + diff_z_ur)./rClass + diff_zz_uz) + mu * ( (rClass .* diff_z_uz + uz_00)./rClass  + diff_zz_uz);

    u0dotdotClass_r = (1/rho) *( forceClass_r  + extforce_r) ;
    u0dotdotClass_z = (1/rho) *( forceClass_z  + extforce_z) ;

    %calculate u0dot
    u0dotClass_r = uolddotClass_r + dt * 0.5 * uolddotdotClass_r + dt * 0.5 * u0dotdotClass_r;
    u0dotClass_z = uolddotClass_z + dt * 0.5 * uolddotdotClass_z + dt * 0.5 * u0dotdotClass_z;
    
    %% For the next loop
    uoldClass_r = u0Class_r;
    uoldClass_z = u0Class_z;

    uolddotClass_r = u0dotClass_r;
    uolddotClass_z = u0dotClass_z;

    uolddotdotClass_r = u0dotdotClass_r;
    uolddotdotClass_z = u0dotdotClass_z;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Bond breaking criteria
    %%Nbd = ((S - snot ) < 0).* Nbd;
    %Nbd = ((stretch - snot ) < 0).* Nbd;
    
    %calculate u0dot
    u0dot = uolddot + dt * 0.5 * uolddotdot + dt * 0.5 * u0dotdot;
    
    uold = u0;
    uolddot = u0dot;
    uolddotdot = u0dotdot;
    
    %loop ends
	fprintf(repmat('\b',1,10));
	fprintf('t=%d',t);

	if mod(t, 50) == 0
		savenewpos(Nbd, u0, 1, dx, dy, nx, ny, imgcounter, f, 'peri_');

		%% save the solution to classical problem
		savenewpos(Nbd, [reshape(u0Class_r', [],1) reshape(u0Class_z', [],1)], 1, dx, dy, nx, ny, imgcounter, f, 'classical_');


		imgcounter = imgcounter + 1;
	end
%     if mod(t,50) == 0
%         drawmesh(Nbd);
%         M(ind) = getframe(f);
%         ind = ind +1;
%     end

	%if t ==100
	%	disp('first layer')
	%	u0(1:5,:)

	%	disp('third layer')
	%	u0(403:408,:)

	%	disp('totalintforce')
	%	1/rho * totalintforce(1:5,:)
	%	1/rho * extforce(1:5,:)

	%	1/rho * totalintforce(403:408,:)
	%	1/rho * extforce(403:408,:)
	%	
	%	%return
	%end

end
toc
out = Nbd;
drawmesh(out, dx, dy, nx, ny);
%movarr = M;



