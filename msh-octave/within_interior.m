function [is_inside] = within_interior(A_min, A_max, B_min, B_max) 
% Determines if the box A lies entirely within the interior of box B
% Input:
%	A_min: 
%	A_max: 
%	: 
%
% Output:
%	out:

A_Left = A_min(1);
A_Right = A_max(1);
A_Top = A_max(2);
A_Bottom = A_min(2);

B_Left = B_min(1);
B_Right = B_max(1);
B_Top = B_max(2);
B_Bottom = B_min(2);

is_inside =  ( B_Left < A_Left && A_Right < B_Right && B_Top > A_Top && B_Bottom < A_Bottom );
