% cross-correlation of vectors with delay r
% caution: A and B must be row vectors of equal length
function out = crossc(A, B, antAPos, antBPos, r)

speed = 1e-2;

D = 5;	% irrelevant since it gets cancelled

% find dtau from r
tauA = 1/speed *( D - dot(antAPos, r));
tauB = 1/speed *( D - dot(antBPos, r));

% making dtau an integer
dtau = floor(tauA - tauB)

out =	sum(A .* circshift(B, dtau, 2));

