function out = shiftMat(M, i, j, nx, ny)
%% Shifts a matrix by -i rows and -j columns but
%% zeros out the offset 
%% Caution: nx is the number of columns
%%	    ny is the number of rows

%% Start with a zero vector
zero_rows = zeros(abs(i),nx);
zero_cols = zeros(ny,abs(j));



out = circshift(M, [-i, -j]);

%% replace by zeros
if i > 0
	out(ny-i+1:ny,:) = zero_rows;
elseif i < 0	% if i==0, do nothing
	out(1:-i,:) = zero_rows; % -i is positive
end

if j > 0
	out(:,nx-j+1:nx) = zero_cols;
elseif j < 0	% if j==0, do nothing
	out(:, 1:-j) = zero_cols;
end









