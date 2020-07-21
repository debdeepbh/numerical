function [total_contact_force, total_friction_force, total_damping_force] = pairwise_forces(CurrPos_center, CurrPos_neighbor, velocity_center, velocity_neighbor, Vol_neighbor, contact_radius, normal_stiffness, friction_coefficient, damping_ratio, rho)

% Computes the pairwise forces: contact force, friction force on the central body due to the neighboring body
% 
% Input:
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


%	contact_NbdArr: matrix containing the indices within the neighboring body that are contact-radius close to the nodes of the central body
% nodes in the neighboring body that are contact_radius-close to the nodes in the central body
%% edit this function to output the difference between the points as well!!
[contact_NbdArr] = gen_NbdArr_varlength(CurrPos_center, CurrPos_neighbor, contact_radius, 0);

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

% relative velocity of the center, i.e. the velocity of the nodes of the central body observed from static neighboring nodes
rel_vc_1 = -(vel_neighbor_1(contact_NbdArr + ~contact_NbdArr) - vel_center_1);
rel_vc_2 = -(vel_neighbor_2(contact_NbdArr + ~contact_NbdArr) - vel_center_2);

rel_vc_norm = sqrt(rel_vc_1.^2 + rel_vc_2.^2);
% unit relative velocity of the center
rel_vc_unit_1 = rel_vc_1 ./ (rel_vc_norm + ~rel_vc_norm) .* mask;
rel_vc_unit_2 = rel_vc_2 ./ (rel_vc_norm + ~rel_vc_norm) .* mask;

% (v.e), projection of unit relative velocity of center onto the direction unit vector (from center to neighbor)
rel_vc_normal_projection = rel_vc_unit_1 .* nbd_direction_unit_1 + rel_vc_unit_2 .* nbd_direction_unit_2;

% v - (v.e)e, relative velocity of center in the tangential direction
rel_vc_tangential_1 = (rel_vc_unit_1 - rel_vc_normal_projection .* nbd_direction_unit_1) .* mask;
rel_vc_tangential_2 = (rel_vc_unit_2 - rel_vc_normal_projection .* nbd_direction_unit_2) .* mask;

% |v - (v.e)e|
rel_vc_tangential_norm = sqrt(rel_vc_tangential_1.^2 + rel_vc_tangential_2.^2);

% normalized (v - (v.e)e), unit vector normal to contact vector
rel_vc_tangential_unit_1 = rel_vc_tangential_1 ./ (rel_vc_tangential_norm + ~rel_vc_tangential_norm) .* mask;
rel_vc_tangential_unit_2 = rel_vc_tangential_2 ./ (rel_vc_tangential_norm + ~rel_vc_tangential_norm) .* mask;

% % tangential friction force acts in the opposite direction of velocity

% % original expression (31) of Peridem.pdf, friction force acts normal to contact dir
% friction_force_1 = -friction_coefficient .* contact_force_norm .* rel_vc_tangential_1;
% friction_force_2 = -friction_coefficient .* contact_force_norm .* rel_vc_tangential_2;

% modification, with normalized unit vector
friction_force_1 = -friction_coefficient .* contact_force_norm .* rel_vc_tangential_unit_1;
friction_force_2 = -friction_coefficient .* contact_force_norm .* rel_vc_tangential_unit_2;


% % modification, direction of friction is opposite to the relative velocity
% % magnitude of friction is the norm of normal contact force
%friction_force_1 = -friction_coefficient .* contact_force_norm .* rel_vc_unit_1; 
%friction_force_2 = -friction_coefficient .* contact_force_norm .* rel_vc_unit_2; 

% % magnitude: projection of contact force perpendicular to relative velocity
% % v_perp = e - (e.v)v
% v_perp_1 = (nbd_direction_unit_1 - rel_vc_normal_projection .* rel_vc_unit_1) .* mask; 
% v_perp_2 = (nbd_direction_unit_2 - rel_vc_normal_projection .* rel_vc_unit_2) .* mask; 
% 
% v_perp_norm = sqrt(v_perp_1.^2 + v_perp_2.^2);
% % unit vector perpendicular to velocity, outward of the central body
% v_perp_unit_1 = v_perp_1 ./ (v_perp_norm + ~v_perp_norm) .* mask;
% v_perp_unit_2 = v_perp_2 ./ (v_perp_norm + ~v_perp_norm) .* mask;
% 
% contact_force_projection = (-nbd_contact_force_1) .* v_perp_unit_1 + (-nbd_contact_force_2) .* v_perp_unit_2;

% friction_force_1 = - friction_coefficient .* contact_force_projection .* rel_vc_unit_1;
% friction_force_2 = - friction_coefficient .* contact_force_projection .* rel_vc_unit_2;

total_friction_force = [sum(friction_force_1, 2), sum(friction_force_2,2)];  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            Damping force                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rel_vc_full_proj = rel_vc_1 .* nbd_direction_unit_1 + rel_vc_2 .* nbd_direction_unit_2;

% % relative velocity of the neighbor
% (v.e)e where v is the velocity of the center
rel_v_neighbor_component_1 = rel_vc_full_proj .* nbd_direction_unit_1;
rel_v_neighbor_component_2 = rel_vc_full_proj .* nbd_direction_unit_2;

damping_coefficient = 2 * damping_ratio .* sqrt(normal_stiffness * rho) .* sqrt(cnbd_Vol_neighbor);

damping_force_1 = - damping_coefficient .* rel_v_neighbor_component_1 .* mask;
damping_force_2 = - damping_coefficient .* rel_v_neighbor_component_2 .* mask;

total_damping_force = [sum(damping_force_1, 2), sum(damping_force_2,2)];  
