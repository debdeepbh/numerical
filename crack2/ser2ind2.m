function j = ser2ind2(y)
global nx;
global ny;


%if ((y <= (nx*ny)) .* (1<=y))
    j = mod(y,nx);

    % if there are zero elements, replace them with nx
    j = j + (j==0)*nx;

    %if j ==0
    %    j = nx;
    %end
%else
%     disp('Bad value');
%end
