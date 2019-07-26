function y = ser2pos2(ser, dy, nx)
% get the y value, given the serial
% yes, both ser2pos1 and se2pos1 need dy and nx
% since nx is enough to determine i and j both by 
% looking at the remainder and result of the division
	i = ser2indi(ser, nx);
	y = indi2pos2(i, dy);
