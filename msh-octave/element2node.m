% Input: a matrix whose row index represents facet index, and each row contains the indices of the vertices of that facet 
% Output: A matrix whose row index represents vertex index, and each row contains the neighboring facet indices
function outArr = element2node(input)

max_nbds = 3;	% maximum possible neighboring facets for each vertex

max_vert_ind = max(max(input));

outArr = zeros(max_vert_ind, max_nbds);

for i = 1:max_vert_ind
	[rows, cols] = find(input == i);
	outArr(i, 1:length(rows)) = rows';	% copy the row numbers into the ith row
end




