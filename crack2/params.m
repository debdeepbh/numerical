function [lam_00, lam_10, lam_20] = params(i,j)
% computes the parameters between the points with serial i and j
global delta


posI = ser2pos(i);
r = posI(1);
z = posI(2);

posJ = ser2pos(j);
r_p = posJ(1);
z_p = posJ(2);

a = r.^2 + r_p.^2+ (z_p - z).^2;
b = (-2).* r.* r_p;
% influence function 
fun_00 = @(phi, a, b) ( 1./ (( a +  b .* cos(phi) ).^(3/2)) );
fun_01 = @(phi, a, b) ( cos(phi)./ ((a + b .* cos(phi) ).^(3/2)) );
fun_20 = @(phi, a ,b) ( (cos(phi).^2)./ ( ( a + b .* cos(phi)).^(3/2)) );

K = (delta^2 - (z_p - z).^2 - r_p.^2 - r.^2)./((-2)*r.*r_p);
alpha = acos(K);
if K > 1
	i
	j
	K
	alpha	
	delta
	return
end


lam_00 = integral(@(phi)fun_00(phi, a, b), (-1)*alpha, alpha);
lam_10 = integral(@(phi)fun_01(phi, a, b), (-1)*alpha, alpha);
lam_20 = integral(@(phi)fun_20(phi, a, b), (-1)*alpha, alpha);
