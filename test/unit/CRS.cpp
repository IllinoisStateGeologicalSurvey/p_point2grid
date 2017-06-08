#include <iostream>
#include <math.h>
#include "CRS.hpp"
#include "Point.hpp"

int main(int argc, char * argv) {
	
	struct point *pt, *pt2;
	pt = new point(0.0, 0.0, 0.0);
	pt2 = new point(0.0, 0.0, 0.0);
	CRS *crs = new CRS();
	pt->randomize();
	pt->print();
	pt2->x = pt->x;
	pt2->y = pt->y;
	pt2->z = pt->z;
	crs->transform(pt);
	pt->print();
	crs->reverse(pt);
	pt->print();
	double distance = pt->distance(pt2);
	std::cout << "Difference: " << distance << "\n" << std::flush;
	delete pt;
	delete pt2;
	delete crs;
	return 0;
} 
