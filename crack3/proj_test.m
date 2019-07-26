function [out, u0] = proj_test(initial, dx, dy, nx, ny, delta)
%% Success

%%% adding these
%global length;
%global width;
%
%global dx;
%global dy;
%global delta;
%
%%% computing these to make the available for ser2pos etc
%global nx;
%global ny;
%
%
%global Nbd;

Nbd = initial;

imgcounter = 1;

%global Xdiff1;
%global Xdiff2;
%global Xdiffnorm;
%

%nx = length/dx+1;
%ny = width/dy+1;
totalnodes = nx*ny;

%ind =1;

%initial condition
sigma=12e6/dx;
rho=2440;
Gnot = 135;
E = 72e9;
nu = 0.22;       % commnet on paper: the effective poisson ration is 1/3

snot = sqrt(4 * pi * Gnot /(9*E*delta));
cnot = 6*E/( pi * (delta^3) * (1 - nu));

f = figure('visible','off');

%internal force matrices
[m, n] = size(Nbd);

Udiff1 = zeros(m,n);
Udiff2 = zeros(m,n);

extforce = zeros(totalnodes,2); % extforce(i,:) = external force on i-th serial
%% Since our y-axis extends downward, the external forces directions are flipped in the $y$ direction.
for j=1:nx
    %extforce(ind2ser(1,j,nx, ny),2) = -(sigma);	%% flipped y direction
    %extforce(ind2ser(ny,j,nx, ny),2) = sigma;
    extforce(ind2ser(1,j,nx, ny),2) = -(sigma);	
    extforce(ind2ser(ny,j,nx, ny),2) = sigma;
end



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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 25e-9;
tic
for t = 1:1800

    % loop starts
    u0 = uold + dt * uolddot + dt * dt * 0.5 * uolddotdot;
    
    restrict = (Nbd > 0);
    
    %for i = 1:totalnodes
    %    for j=1:n
    %        if Nbd(i,j) > 0
    %            Udiff1(i,j) = u0(Nbd(i,j),1) - u0(i,1);
    %            Udiff2(i,j) = u0(Nbd(i,j),2) - u0(i,2);
    %        end
    %    end
    %end
    %
    %Bigvect1 = Udiff1 + Xdiff1;
    %Bigvect2 = Udiff2 + Xdiff2;
    %
    %Bigvectnorm = sqrt(Bigvect1.*Bigvect1 + Bigvect2.*Bigvect2);    % norm(ucap - u + xacp -x)
    %S = (Bigvectnorm - Xdiffnorm ) ./ Xdiffnorm;    %gives NaN when an element  is not a nbd
    
    %Multipl = cnot * (S ./ Bigvectnorm);
    %Force1 = Multipl .* Bigvect1;
    %Force2 = Multipl .* Bigvect2;

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

    eta_p_xi_r = eta_r + xi_r;
    eta_p_xi_z = eta_z + xi_z;

    eta_p_xi_norm = sqrt(eta_p_xi_r.^2 + eta_p_xi_z.^2);

    % to avoid dividing by zero 
	stretch = (eta_p_xi_norm - xi_norm)./(xi_norm + (~xi_norm));

    unitvec_r =  eta_p_xi_r ./ eta_p_xi_norm;
    unitvec_z =  eta_p_xi_z ./ eta_p_xi_norm;

    Force1 = cnot.* stretch .* unitvec_r;
    Force2 = cnot.* stretch .* unitvec_z;

    %% mask it
    Force1 = Force1 .* (~~Nbd);
    Force2 = Force2 .* (~~Nbd);

    
    totalintforce = [sum(Force1.*restrict, 2) sum(Force2.*restrict,2)];    %sum(A) ad totalintforce2 = sum(Force2, 2);
    
    u0dotdot = (1/rho) *( totalintforce * (dx*dy) + extforce) ;
    %u0dotdot = (1/rho) *( totalintforce * (dx*dy)) + extforce ;
    
    
    %% Bond breaking criteria
    %%Nbd = ((S - snot ) < 0).* Nbd;
    %Nbd = ((stretch - snot ) < 0).* Nbd;
    
    %calculate u0dot
    u0dot = uolddot + dt * 0.5 * uolddotdot + dt * 0.5 * u0dotdot;
    
    uold = u0;
    uolddot = u0dot;
    uolddotdot = u0dotdot;
    
    %loop ends
    t
%     if mod(t,50) == 0
%         drawmesh(Nbd);
%         M(ind) = getframe(f);
%         ind = ind +1;
%     end
	if mod(t, 50) == 0
		savenewpos(Nbd, u0, 2, dx, dy, nx, ny, imgcounter, f, 'test');
		imgcounter = imgcounter + 1;
	end
end
toc
out = Nbd;
drawmesh(out, dx, dy, nx, ny);
%movarr = M;



