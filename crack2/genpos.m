function genpos(Nbd)
% generate position matrices in the reference config
% following the same rule as the neighborhood matrix
global Pos1
global Pos2

% these vectors of length `totalnodes` store the position of all nodes
global PosVec1
global PosVec2

ss = size(Nbd);

PosVec1 = ser2pos1( (1:ss(1))'); % want it to be a column vector, thus the '
PosVec2 = ser2pos2( (1:ss(1))'); % want it to be a column vector, thus the '


missing1 = ~Nbd;	% this is the matrix with 1's in the missing positions
			% hence ~missing1 is the matrix with 0 in the missing positions

% This line populates the missing serials by 1, evaluates the Position1, then sets the missing positions to zero
Pos1 = ser2pos1(Nbd + missing1).* (~missing1);
Pos2 = ser2pos2(Nbd + missing1).* (~missing1);
