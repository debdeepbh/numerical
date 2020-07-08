function [total_contact_force, nbd_contact_force_1, nbd_contact_force_2, nbd_direction_unit_1, nbd_direction_unit_2] = contact_force(contact_NbdArr, CurrPos_center, CurrPos_neighbor, Vol_neighbor, contact_radius, normal_stiffness) 
% Computes the contact force on the central body due to the neighboring body
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
