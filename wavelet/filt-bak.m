% case statement not working, syntax error
% filters
function [u, v] = filt(TYPE, N)
switch TYPE
case 'shan'
	% shannon filters
	if mod(N,4) != 0
		printf('error:dimension not divisible by 4');
	else
		ucap = [sqrt(2).*ones(1,N/4), i, zeros(1,N/2-1), -i, sqrt(2)*ones(1,N/4 -1)];
		vcap = [zeros(1,N/4), 1, sqrt(2).*ones(1,N/2-1), 1, zeros(1,N/4 -1)];
		u = real(ifft(ucap));
		v = real(ifft(vcap));
		%v = getother(u);
	end
case { 'd2', 'haar' }
	u = [1/sqrt(2), 1/sqrt(2), zeros(1,N-2)];
	v = [1/sqrt(2), -1/sqrt(2), zeros(1,N-2)];
case 'd4'
	a = 1 + sqrt(3);
	b = 3 + sqrt(3);
	c = 3 - sqrt(3);
	d = 1 - sqrt(3);
	mult = sqrt(2)/8;
	u = mult * [a, b, c, d, zeros(1,N-4)];
	v = mult * [-b, a, zeros(1,N-4), -d, c];
case 'd6'
	% Daubechies filters
	a = 1 - sqrt(10);
	b = 1 + sqrt(10);
	c = sqrt(5+2*sqrt(10));

	% without the parenthesis, dimension of u increased by 1!
	u = (sqrt(2)/32).*[ (b+c), (2*a + 3*b + 3*c), (6*a + 4*b + 2*c), (6*a + 4*b -2*c), (2*a + 3*b-3*c), (b-c), zeros(1,N-6)];
	% adding 1 to the indices in definition 
	v = [ -u(2), u(1), zeros(1,N-6), -u(6), u(5), -u(4), u(3) ];
case 'd8'
	% not normalized u
	u =[0.32580343, 1.01094572, 0.89220014, -0.03957503, -0.26450717, 0.0436163, 0.0465036, -0.01498699, zeros(1,N-8)];
	v = getother(u);
	% normalize the vectors
	u = u/norm(u);
	v = u/norm(v);

case 'd10'
	u = [ 0.22641898, 0.85394354, 1.02432694, 0.19576696, -0.34265671, -0.04560113, 0.10970265, -0.0088268, -0.01779187, 0.00471742793, zeros(1,N-10)];

	v = getother(u);
	% normalize the vectors
	u = u/norm(u);
	v = u/norm(v);
case 'd12'

	u = [0.15774243, 0.69950381, 1.06226376, 0.44583132, -0.3199866, -0.18351806, 0.13788809, 0.03892321, -0.04466375, 0.000783251152, 0.00675606236, -0.00152353381, zeros(1,N-12)];
	v = getother(u);
	% normalize the vectors
	u = u/norm(u);
	v = u/norm(v);

case 'd14'

	u = [0.11009943, 0.56079128, 1.03114849, 0.66437248, -0.20351382, -0.31683501, 0.1008467, 0.11400345, -0.05378245, -0.02343994, 0.01774979, 0.000607514995, -0.00254790472, 0.000500226853, zeros(1,N-14)];
	v = getother(u);
	% normalize the vectors
	u = u/norm(u);
	v = u/norm(v);

case 'd16'

	u = [0.07695562, 0.44246725, 0.95548615, 0.82781653, -0.02238574, -0.40165863, 0.000668194092, 0.18207636, -0.0245639, -0.06235021, 0.01977216, 0.01236884, -0.00688771926, -0.000554004549, 0.000955229711, -0.000166137261, zeros(1,N-16)];
	v = getother(u);
	% normalize the vectors
	u = u/norm(u);
	v = u/norm(v);

case 'd18'

	u = [0.05385035, 0.3448343, 0.85534906, 0.92954571, 0.18836955, -0.41475176, -0.13695355, 0.21006834, 0.043452675, -0.09564726, 0.000354892813, 0.03162417, -0.00667962023, -0.00605496058, 0.00261296728, 0.000325814671, -0.000356329759, 5.5645514e-05, zeros(1,N-18)];
	v = getother(u);
	% normalize the vectors
	u = u/norm(u);
	v = u/norm(v);

case 'd20'

	u = [ 0.03771716, 0.26612218, 0.74557507, 0.97362811, 0.39763774, -0.3533362, -0.27710988, 0.18012745, 0.13160299, -0.10096657, -0.04165925, 0.04696981, 0.00510043697, -0.015179, 0.00197332536, 0.00281768659, -0.00096994784, -0.000164709006, 0.000132354367, -1.875841e-05, zeros(1,N-20)];
	v = getother(u);
	% normalize the vectors
	u = u/norm(u);
	v = u/norm(v);
otherwise
	printf('Type of filter not valid.\n');
end
