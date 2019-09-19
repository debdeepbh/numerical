function j = ser2indj(ser,nx)

%size(ser)
%size(nx)
    j = mod(ser,nx);

    %if j ==0
    %    j = nx;
    %end
    % replace all the zero values by nx
    j = j + (~j)*nx;
    
