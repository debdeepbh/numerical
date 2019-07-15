function x = ser2pos1(serial)
% outputs the x-coordinate
global dx;


	j = ser2ind2(serial);

	% j is the column number, associated with the x-value
	x = (j-1)*dx + 0.1;
	%x = (j-1)*dx;

