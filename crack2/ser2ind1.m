function i = ser2ind1(y)

% outputs the row number, given the serial

global nx;
global ny;

%if ((y <= (nx*ny)) .* (1<=y))
	i = floor((y-1)/nx)+1;
%else
%     disp('Bad value');
%end
