lcCircle = 50;
lcSquare = 150;
lcAvg    = (lcCircle+lcSquare)/2;
//lc = 50;

y0 = 0;
y1 = 2000;
x0 = 0;
x1 = 2000;


//square
Point(1) = {x0, y0, 0, lcSquare};
Point(2) = {x1, y0, 0, lcSquare};
Point(3) = {x1, y1, 0, lcSquare};
Point(4) = {x0, y1, 0, lcSquare};
Line(1) = {1,2};
Line(2) = {3,2};
Line(3) = {3,4};
Line(4) = {4,1};

Periodic Line {2} = {-4} ;

Line Loop(5) = {1,-2,3,4};


//circle
ccx = 1000;
ccy = 1000;

rad = 250;

c1x = ccx;
c1y = ccy+rad; 
c2x = ccx+rad; 
c2y = ccy; 
c3x = ccx;
c3y = ccy-rad;
c4x = ccx-rad;
c4y = ccy;

Point(10) = {ccx,ccy, 0, lcCircle};
Point(11) = {c1x,c1y, 0, lcCircle};
Point(12) = {c2x,c2y, 0, lcCircle};
Point(13) = {c3x,c3y, 0, lcCircle};
Point(14) = {c4x,c4y, 0, lcCircle};

Circle(21) = {11,10,12};
Circle(22) = {12,10,13};
Circle(23) = {13,10,14};
Circle(24) = {14,10,11};

Line Loop(30) = {-21,-22,-23,-24};


Plane Surface(1) = {5, 30} ; //EXTERIOR
Plane Surface(2) = {30};     //INTERIOR


Field[1] = Attractor;
Field[1].NodesList = {10};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lcCircle;
Field[2].LcMax = lcSquare;
Field[2].DistMin = rad;
Field[2].DistMax = rad*(1+1/2);

Background Field = 2;

// Don't extend the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 0;

