function y = ser2pos2(serial)
% outputs the y-coordinate
global dy;

	i = ser2ind1(serial);

	% i is the row number, associated with the y-value
	y = (i-1)*dy;

