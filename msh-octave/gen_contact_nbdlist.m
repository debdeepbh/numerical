function [contact_NbdArr] = gen_contact_nbdlist(CurrPos_multi, contact_radius) 
% Returns the node neighborhood array associated with neighboring particles
% n = total nodes in a particle
% total_particles = total number of particles
% Input:
%	CurrPos_multi: Current position of nodes of all the particles, nx2xtotal_particles
%	contact_radius: contact radius, scalar
%
% Output:
%	contact_NbdArr: neighborhood array that contains the indices of nodes inside another particle that are within the contact radius, n x l x choose(total_particles, 2), l = max number of nodes within contact radius


[total_nodes, dimension, total_particles]  = size(CurrPos_multi);
indices_all = (1:total_nodes)';

% generate contact neighbor list
contact_indices = nchoosek(1:total_particles, 2);
[total_contacts, temp] = size(contact_indices);



%contact_NbdArr = zeros(total_contact, 


nbd_count = zeros(total_nodes,1);


for k = 1:total_contacts

    [particle_center, particle_neighbor] = contact_indices(k,:);

    for i = 1:total_nodes	% i is a node in particle_center
	diff = CurrPos(:,:,particle_neighbor) - CurrPos(i,:, particle_center);
	neighbor_indices = ( sqrt(diff(:,1).^2 + diff(:,2).^2)  < contact_radius) .* indices_all;

	nbdlist{i} = nonzeros(neighbor_indices);
	nbd_count(i) = length(nbdlist{i});
    end

    % Alternative
    diffMat = repmat(CurrPos(:,:,particle_neighbor), 1, total_nodes) - reshape(CurrPos(:,:,particle_center)', 1, []);
    neighbor_indices = ( sqrt(diffMat(:, 1:2:end).^2 + diffMat(:, 2:2:end).^2) < contact_radius ) .* indices_all;


    % convert to matrix from cell array
    NbdArr = zeros(total_nodes, max(nbd_count));
    for i = 1:total_nodes
	    NbdArr(i, 1:nbd_count(i)) = nbdlist{i};
    end

end




