% case statement not working, syntax error
% filters
function [u, v] = filt(TYPE, N)
switch TYPE

case 'bior13d'	% decomposition bi-orthogoanl 1.3
	% these values are taken from python library
	% http://wavelets.pybytes.com/wavelet/bior1.3/
	% the adjustment is done to make it fit my notations
	% here is a recipe to adapt the data from python library to
	% my notations
	u = [-0.08838834764831845 0.08838834764831845 0.7071067811865476 0.7071067811865476 0.08838834764831845 -0.08838834764831845];
	v = [0.0 0.0 -0.7071067811865476 0.7071067811865476 0.0 0.0];
	%u = shift([zeros(1,N-6) u],1);
	%v = shift([v zeros(1,N-6)],-1);

case 'bior13r'	% reconstruction bi-orthogonal 1.3
	u = [0.0 0.0 0.7071067811865476 0.7071067811865476 0.0 0.0];
	v = [-0.08838834764831845 -0.08838834764831845 0.7071067811865476 -0.7071067811865476 0.08838834764831845 0.08838834764831845];
	%u = shift([u zeros(1,N-6)],0);
	%v = shift([zeros(1,N-6) v],2);
	%u = ifft(conj(fft(u)));
	%v = ifft(conj(fft(v)));
case 'bior33d'
	u = [0.06629126073623884 -0.19887378220871652 -0.15467960838455727 0.9943689110435825 0.9943689110435825 -0.15467960838455727 -0.19887378220871652 0.06629126073623884];
	v = [0.0 0.0 -0.1767766952966369 0.5303300858899107 -0.5303300858899107 0.1767766952966369 0.0 0.0];
case 'bior33r'
	u = [0.0 0.0 0.1767766952966369 0.5303300858899107 0.5303300858899107 0.1767766952966369 0.0 0.0];
	v = [0.06629126073623884 0.19887378220871652 -0.15467960838455727 -0.9943689110435825 0.9943689110435825 0.15467960838455727 -0.19887378220871652 -0.06629126073623884];


case 'testd6d'	% from python for testing purposes
	u = [0.035226291882100656 -0.08544127388224149 -0.13501102001039084 0.4598775021193313 0.8068915093133388 0.3326705529509569];
        v = [-0.3326705529509569 0.8068915093133388 -0.4598775021193313 -0.13501102001039084 0.08544127388224149 0.035226291882100656]; 
	u = shift([zeros(1,N-6) u],1);
	v = shift([v zeros(1,N-6)],-1);
case 'testd6r'
	u = [0.3326705529509569 0.8068915093133388 0.4598775021193313 -0.13501102001039084 -0.08544127388224149 0.035226291882100656]; 
	v = [0.035226291882100656 0.08544127388224149 -0.13501102001039084 -0.4598775021193313 0.8068915093133388 -0.3326705529509569]; 

	u = shift([u zeros(1,N-6)],0);
	v = shift([zeros(1,N-6) v],2);
	u = ifft(conj(fft(u)));
	v = ifft(conj(fft(v)));


case 'meyer'
	u = [-1.009999956941423e-12 8.519459636796214e-09 -1.111944952595278e-08 -1.0798819539621958e-08 6.066975741351135e-08 -1.0866516536735883e-07 8.200680650386481e-08 1.1783004497663934e-07 -5.506340565252278e-07 1.1307947017916706e-06 -1.489549216497156e-06 7.367572885903746e-07 3.20544191334478e-06 -1.6312699734552807e-05 6.554305930575149e-05 -0.0006011502343516092 -0.002704672124643725 0.002202534100911002 0.006045814097323304 -0.006387718318497156 -0.011061496392513451 0.015270015130934803 0.017423434103729693 -0.03213079399021176 -0.024348745906078023 0.0637390243228016 0.030655091960824263 -0.13284520043622938 -0.035087555656258346 0.44459300275757724 0.7445855923188063 0.44459300275757724 -0.035087555656258346 -0.13284520043622938 0.030655091960824263 0.0637390243228016 -0.024348745906078023 -0.03213079399021176 0.017423434103729693 0.015270015130934803 -0.011061496392513451 -0.006387718318497156 0.006045814097323304 0.002202534100911002 -0.002704672124643725 -0.0006011502343516092 6.554305930575149e-05 -1.6312699734552807e-05 3.20544191334478e-06 7.367572885903746e-07 -1.489549216497156e-06 1.1307947017916706e-06 -5.506340565252278e-07 1.1783004497663934e-07 8.200680650386481e-08 -1.0866516536735883e-07 6.066975741351135e-08 -1.0798819539621958e-08 -1.111944952595278e-08 8.519459636796214e-09 -1.009999956941423e-12 0.0];
	u = [u zeros(1,N-length(u))];
	v = getother(u);

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

case 'd10'
	u = [ 0.22641898, 0.85394354, 1.02432694, 0.19576696, -0.34265671, -0.04560113, 0.10970265, -0.0088268, -0.01779187, 0.00471742793, zeros(1,N-10)];

case 'd12'

	u = [0.15774243, 0.69950381, 1.06226376, 0.44583132, -0.3199866, -0.18351806, 0.13788809, 0.03892321, -0.04466375, 0.000783251152, 0.00675606236, -0.00152353381, zeros(1,N-12)];

case 'd14'

	u = [0.11009943, 0.56079128, 1.03114849, 0.66437248, -0.20351382, -0.31683501, 0.1008467, 0.11400345, -0.05378245, -0.02343994, 0.01774979, 0.000607514995, -0.00254790472, 0.000500226853, zeros(1,N-14)];

case 'd16'

	u = [0.07695562, 0.44246725, 0.95548615, 0.82781653, -0.02238574, -0.40165863, 0.000668194092, 0.18207636, -0.0245639, -0.06235021, 0.01977216, 0.01236884, -0.00688771926, -0.000554004549, 0.000955229711, -0.000166137261, zeros(1,N-16)];

case 'd18'

	u = [0.05385035, 0.3448343, 0.85534906, 0.92954571, 0.18836955, -0.41475176, -0.13695355, 0.21006834, 0.043452675, -0.09564726, 0.000354892813, 0.03162417, -0.00667962023, -0.00605496058, 0.00261296728, 0.000325814671, -0.000356329759, 5.5645514e-05, zeros(1,N-18)];

case 'd20'

	u = [ 0.03771716, 0.26612218, 0.74557507, 0.97362811, 0.39763774, -0.3533362, -0.27710988, 0.18012745, 0.13160299, -0.10096657, -0.04165925, 0.04696981, 0.00510043697, -0.015179, 0.00197332536, 0.00281768659, -0.00096994784, -0.000164709006, 0.000132354367, -1.875841e-05, zeros(1,N-20)];
otherwise
	printf('Type of filter not valid.\n');
end

% normalizing the vectors, taken from wikipedia
switch TYPE
	case { 'd8', 'd10', 'd12', 'd14', 'd16', 'd18', 'd20' }
		v = getother(u);
		u = u/sqrt(2);
		v = v/sqrt(2);
end

% modification to be done when copying from python library blindly
switch TYPE	% deconstruction vectors
case { 'bior13d', 'bior33d'}
	u = shift([zeros(1,N-length(u)) u],1);
	v = shift([v zeros(1,N-length(v))],-1);
end
switch TYPE	% reconstruction vectors
case { 'bior13r', 'bior33r'}
	u = shift([u zeros(1,N-length(u))],0);
	v = shift([zeros(1,N-length(v)) v],2);

	% Find the tilde of that
	u = ifft(conj(fft(u)));
	v = ifft(conj(fft(v)));
end



