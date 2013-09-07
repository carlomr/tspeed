lcFine   = 50;
lcCoarse = 300;
lcAvg    = (lcFine+lcCoarse)/2;
//lc = 50;

y0 = 0;
y1 = 2000;
x0 = 0;
x1 = 2000;


Point(1) = {x0, y0, 0, lcFine};
Point(2) = {x1, y0, 0, lcCoarse} ;
Point(3) = {x1, y1, 0, lcFine} ;
Point(4) = {x0, y1, 0, lcFine} ;

Line(1) = {1,2} ;
Line(2) = {3,2} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;
Line(5) = {1,3} ;

Line Loop(5) = {4,5,3} ;
Line Loop(6) = {1,-2,-5};
Plane Surface(2) = {5} ;
Plane Surface(1) = {6};

