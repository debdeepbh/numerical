function [NbdArr_out, u0] = simulateMultiple(Pos, NbdArr, Vol, extforce, delta, xi_1, xi_2, xi_norm)


% break bonds or not
break_bonds = 0;

close all

% get the material properties
[delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta);


total_nodes = length(NbdArr);


% Initial data
uold = zeros(total_nodes,2);
uolddot = zeros(total_nodes,2);
uolddotdot = zeros(total_nodes,2);

longvec = (1:total_nodes)';

imgcounter = 1;
%dt = 25e-9;
dt = 25e-8;

f = figure('visible','off');
%f = figure('visible','on');
tic

for t = 1:1800
%for t = 1:2800

    % loop starts
    u0 = uold + dt * uolddot + dt * dt * 0.5 * uolddotdot;
    
    restrict = (NbdArr > 0);

    disp_1 = u0(:,1);
    disp_2 = u0(:,2);

    %% Displacements of the centers
    u_1 = disp_1(longvec);
    u_2 = disp_2(longvec);

    %% Displacements of the neighbors
    u_1_p = disp_1(NbdArr + (~NbdArr));
    u_2_p = disp_2(NbdArr + (~NbdArr));

    eta_1 = u_1_p - u_1;
    eta_2 = u_2_p - u_2;

    eta_p_xi_1 = eta_1 + xi_1;
    eta_p_xi_2 = eta_2 + xi_2;

    eta_p_xi_norm = sqrt(eta_p_xi_1.^2 + eta_p_xi_2.^2);

    % to avoid dividing by zero 
	stretch = (eta_p_xi_norm - xi_norm)./(xi_norm + (~xi_norm));

    unitvec_1 =  eta_p_xi_1 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* restrict;
    unitvec_2 =  eta_p_xi_2 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* restrict;

    %unitvec_1 =  eta_p_xi_1 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* (~~eta_p_xi_norm);
    %unitvec_2 =  eta_p_xi_2 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* (~~eta_p_xi_norm);

    Force1 = cnot.* stretch .* unitvec_1;
    Force2 = cnot.* stretch .* unitvec_2;

    %% mask it
    Force1 = Force1 .* (~~NbdArr);
    Force2 = Force2 .* (~~NbdArr);


    totalintforce = [sum(Force1.*restrict, 2) sum(Force2.*restrict,2)];    %sum(A) ad totalintforce2 = sum(Force2, 2);


    
    dx = 0.002/4;
    u0dotdot = (1/rho) *( totalintforce .* Vol + extforce) ;
    
%    nbds = nonzeros(NbdArr(1,:));
%    vec = extforce;
%	tif1 = vec(:,1);
%	tif2 = vec(:,2);
%    [tif1(nbds) tif2(nbds)]
%    t

    
    if break_bonds == 1
	    % Bond breaking criteria
	    %NbdArr = ((S - snot ) < 0).* NbdArr;
	    NbdArr = ((stretch - snot ) < 0).* NbdArr;
    end


    
    %calculate u0dot
    u0dot = uolddot + dt * 0.5 * uolddotdot + dt * 0.5 * u0dotdot;

    
    uold = u0;
    uolddot = u0dot;
    uolddotdot = u0dotdot;


    
    %loop ends
%     if mod(t,50) == 0
%         drawmesh(NbdArr);
%         M(ind) = getframe(f);
%         ind = ind +1;
%     end
 % if t == 2
 %       savenewpos2(u0, Pos, 3, imgcounter, f, 'test')
 %       pause
 % end
	if mod(t, 50) == 0
		t
		savenewpos2(u0, Pos, 3, imgcounter, f, 'test')
		imgcounter = imgcounter + 1;
	end

%scaling = 10;
%newPos = u0 + scaling * Pos;
%defnorm = sqrt(u0(:,1).^2 + u0(:,2).^2);
%scatter(newPos(:,1), newPos(:,2), 5, defnorm);
%pause

end

toc

% output
NbdArr_out = NbdArr;

%drawmesh(out, dx, dy, nx, ny);



