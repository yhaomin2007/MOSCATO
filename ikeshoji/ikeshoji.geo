D=1.0; //1E-3;
Point(1) = {0,0,0};
Point(2) = {9*D,0,0};
Point(3) = {9*D,D,0};
Point(4) = {0,D,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1};

Transfinite Curve{1,3} = 81 ;
Transfinite Curve{2,4} = 21 Using Bump 0.1;

Physical Line("bot",1) = {1};
Physical Line("right",2) = {2};
Physical Line("top",3) = {3};
Physical Line("left",4) = {4};

Physical Surface("Fluid") = {1};
