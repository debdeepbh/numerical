 // File created with Octave
// Points
Point(1) = {-0.001,0,0,0.0002};
Point(2) = {0.001,0,0,0.0002};
Point(3) = {0,0.001732050807568877,0,0.0002};
// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};
Line Loop(4) = {1,2,3};
Plane Surface(5) = {4};
