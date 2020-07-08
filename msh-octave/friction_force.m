function [total_friction_force] = friction_force(contact_NbdArr, vel_center_1, vel_center_2, vel_neighbor_1, vel_neighbor_2, nbd_contact_force_1, nbd_contact_force_2, nbd_direction_unit_1, nbd_direction_unit_2, Vol_center, friction_coefficient) 
% Function description
% Input:
%	contact_NbdArr: 
%	vel_center_1: 
%	vel_center_2: 
%	vel_neighbor_1: 
%	vel_neighbor_2: 
%	nbd_contact_force_1: 
%	nbd_contact_force_2: 
%	nbd_direction_unit_1: 
%	nbd_contact_force_2: 
%	Vol_center: 
%	friction_coefficient: 
%
% Output:
%	out:


mask = (~~contact_NbdArr);

% getting the contact force on i-th node from the force density by multiplying by the volume of the node
contact_force_norm = sqrt(nbd_contact_force_1.^2 + nbd_contact_force_2.^2) .* Vol_center;

% relative velocity
vel_direction_1 = (vel_neighbor_1(contact_NbdArr + ~contact_NbdArr) - vel_center_1) .* mask;
vel_direction_2 = (vel_neighbor_2(contact_NbdArr + ~contact_NbdArr) - vel_center_2) .* mask;

vel_direction_norm = sqrt(vel_direction_1.^2  + vel_direction_2.^2);

%  in the direction from center to neighbor
vel_direction_unit_1 = vel_direction_1 ./ (vel_direction_norm + ~vel_direction_norm) .* mask;
vel_direction_unit_2 = vel_direction_2 ./ (vel_direction_norm + ~vel_direction_norm) .* mask;

% projection of unit velocity onto the direction unit vector
vel_normal_projection = vel_direction_unit_1 .* nbd_direction_unit_1 + vel_direction_unit_2 .* nbd_direction_unit_2; 

% v - (v.e)e
vel_tangential_1 = (vel_direction_1 - vel_normal_projection .* nbd_direction_unit_1) .* mask; 
vel_tangential_2 = (vel_direction_2 - vel_normal_projection .* nbd_direction_unit_2) .* mask;

% tangential friction force acts in the opposite direction of velocity
friction_force_1 = friction_coefficient .* contact_force_norm .* vel_tangential_1;
friction_force_2 = friction_coefficient .* contact_force_norm .* vel_tangential_2;

total_friction_force = [sum(friction_force_1, 2), sum(friction_force_2,2)];  


