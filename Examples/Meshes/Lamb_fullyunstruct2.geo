//lcFine   = 50;
//lcCoarse = 300;
lc = 50;

y0 = 0;
y1 = 2000;
x0 = 0;
x1 = 4000;

coarseParam = 8;


Point(1) = {x0, y0, 0, lc*coarseParam};
Point(2) = {x1, y0, 0, lc*coarseParam} ;
Point(3) = {x1, y1, 0, lc*coarseParam} ;
Point(4) = {x0, y1, 0, lc*coarseParam} ;

Line(1) = {1,2} ;
Line(2) = {3,2} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Line Loop(5) = {4,1,-2,3} ;
Plane Surface(6) = {5} ;

//Attractor on top edge
Field[1] = Attractor;
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = {3};

//
// LcMax -                         /------------------
//                               /
//                             /
//                           /
// LcMin -o----------------/
//        |                |       |
//     Attractor       DistMin   DistMax
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc;
Field[2].LcMax = lc*coarseParam;
Field[2].DistMin = (y1-y0)/4;
Field[2].DistMax = (y1-y0)*1/2;

Background Field = 2;


// Don't extend the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 0;



