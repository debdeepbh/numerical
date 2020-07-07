function [does_intersect] = intersects_box(A_min, A_max, B_min, B_max) 
% Check if two rectangles intersect
% Input:
%	min_A: 
%	max_A: 
%	min_B: 
%	max_B: 
%
% Output:
%	does_intersect:

A_Left = A_min(1);
A_Right = A_max(1);
A_Top = A_max(2);
A_Bottom = A_min(2);

B_Left = B_min(1);
B_Right = B_max(1);
B_Top = B_max(2);
B_Bottom = B_min(2);

does_intersect =  (A_Left <= B_Right && A_Right >= B_Left && A_Top >= B_Bottom && A_Bottom <= B_Top );

