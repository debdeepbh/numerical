function x = ser2pos1(ser, dx, nx)
% get the x- value, given the serial
	j = ser2indj(ser, nx);
	x = indj2pos1(j, dx);
