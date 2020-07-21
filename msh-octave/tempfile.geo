 // File created with Octave
// Points
Point(1) = {-0.002,-0.00025,0,0.0002};
Point(2) = {0.002,-0.00025,0,0.0002};
Point(3) = {0.002,0,0,0.0002};
Point(4) = {-0.002,0,0,0.0002};
// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
