function [total_contact_force, total_friction_force] = pairwise_forces(contact_NbdArr, CurrPos_center, CurrPos_neighbor, velocity_center, velocity_neighbor, Vol_neighbor, contact_radius, normal_stiffness, friction_coefficient)

% Computes the pairwise forces: contact force, friction force on the central body due to the neighboring body
% 
% Input:
%	contact_NbdArr: matrix containing the indices within the neighboring body that are contact-radius close to the nodes of the central body
%	CurrPos_center:  Current position vector of central body, nx2
%	CurrPos_neighbor:Current position vector of neighboring body, mx2
%	Vol_neighbor: Volume associated with the nodes of the neighboring body, mx1
%	contact_radius: scalar
%	normal_stiffness: scalar
%
% Output:
%	total_contact_force: Total contact force being applied on the nodes of the central body due to the neighboring body, nx2
%	direction_unit_1: first component of direction vector from the center to neighbor, nxl, l=#neighbors
%	direction_unit_2: second component of direction vector from the center to neighbor, nxl, l=#neighbors

CurrPos_center_1 = CurrPos_center(:,1);
CurrPos_center_2 = CurrPos_center(:,2);

CurrPos_neighbor_1 = CurrPos_neighbor(:,1);
CurrPos_neighbor_2 = CurrPos_neighbor(:,2);

% mask for contact neighbors
mask = ~~contact_NbdArr;

% volume of the contact neighbors
cnbd_Vol_neighbor = Vol_neighbor(contact_NbdArr + ~contact_NbdArr) .* mask;

% direction vectors from the nodes in the central body to the delta-close nodes in the neighboring body
direction_1 = (CurrPos_neighbor_1(contact_NbdArr + ~contact_NbdArr) - CurrPos_center_1) .* mask;
direction_2 = (CurrPos_neighbor_2(contact_NbdArr + ~contact_NbdArr) - CurrPos_center_2) .* mask;

direction_norm = sqrt(direction_1.^2  + direction_2.^2);

% unit vector from center to neighbor
nbd_direction_unit_1 = direction_1 ./ (direction_norm + ~direction_norm) .* mask;
nbd_direction_unit_2 = direction_2 ./ (direction_norm + ~direction_norm) .* mask;

% positive contact radius
cont_rad_contrib = (contact_radius - direction_norm);
positive_cont_rad_contrib = cont_rad_contrib .* (cont_rad_contrib > 0);

%% Check this formula to see if the volume is properly accounted for
% multiplied by the volumes of contact neighbors
% in the negative direction, i.e. force acts on the center in the direction from the neighbors to the center
nbd_contact_force_1 = - normal_stiffness .* positive_cont_rad_contrib .* cnbd_Vol_neighbor .* nbd_direction_unit_1;
nbd_contact_force_2 = - normal_stiffness .* positive_cont_rad_contrib .* cnbd_Vol_neighbor .* nbd_direction_unit_2;

total_contact_force = [sum(nbd_contact_force_1, 2) sum(nbd_contact_force_2,2)];  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     Friction force computation                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vel_center_1 = velocity_center(:,1);
vel_center_2 = velocity_center(:,2);

vel_neighbor_1 = velocity_neighbor(:,1);
vel_neighbor_2 = velocity_neighbor(:,2);

% this is actually the contact force density
contact_force_norm = sqrt(nbd_contact_force_1.^2 + nbd_contact_force_2.^2);

% relative velocity
vel_direction_1 = (vel_neighbor_1(contact_NbdArr + ~contact_NbdArr) - vel_center_1);
vel_direction_2 = (vel_neighbor_2(contact_NbdArr + ~contact_NbdArr) - vel_center_2);
%vel_direction_1 = (vel_neighbor_1(contact_NbdArr + ~contact_NbdArr) - vel_center_1) .* mask;
%vel_direction_2 = (vel_neighbor_2(contact_NbdArr + ~contact_NbdArr) - vel_center_2) .* mask;

vel_direction_norm = sqrt(vel_direction_1.^2  + vel_direction_2.^2) .* mask;

%  in the direction from center to neighbor
vel_direction_unit_1 = vel_direction_1 ./ (vel_direction_norm + ~vel_direction_norm) .* mask;
vel_direction_unit_2 = vel_direction_2 ./ (vel_direction_norm + ~vel_direction_norm) .* mask;

% projection of unit velocity onto the direction unit vector
vel_normal_projection = vel_direction_unit_1 .* nbd_direction_unit_1 + vel_direction_unit_2 .* nbd_direction_unit_2; 

% v - (v.e)e
vel_tangential_1 = (vel_direction_1 - vel_normal_projection .* nbd_direction_unit_1) .* mask; 
vel_tangential_2 = (vel_direction_2 - vel_normal_projection .* nbd_direction_unit_2) .* mask;

vel_tangential_norm = sqrt(vel_tangential_1.^2 + vel_tangential_2.^2);

vel_tangential_unit_1 = vel_tangential_1 ./ (vel_tangential_norm + ~vel_tangential_norm) .* mask;
vel_tangential_unit_2 = vel_tangential_2 ./ (vel_tangential_norm + ~vel_tangential_norm) .* mask;

% % tangential friction force acts in the opposite direction of velocity
friction_force_1 = friction_coefficient .* contact_force_norm .* vel_tangential_unit_1;
friction_force_2 = friction_coefficient .* contact_force_norm .* vel_tangential_unit_2;
%friction_force_1 = friction_coefficient .* contact_force_norm .* vel_tangential_1;
%friction_force_2 = friction_coefficient .* contact_force_norm .* vel_tangential_2;

total_friction_force = [sum(friction_force_1, 2), sum(friction_force_2,2)];  


