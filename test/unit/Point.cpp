#include <iostream>
#include <vector>
#include "QuadTree.hpp"
#include "CRS.hpp"

using namespace std;

int main(int argc, char * argv[]) {

	struct point *pt = new point(10.0, 20.0, 30.0);
	struct point *pt2 = new point(5.0, 15.0, 25.0);
	pt->print();
	pt2->print();
	double dist = pt->distance(pt2);
	double alt = pt->altitd(pt2);
	double azim = pt->azimth(pt2);
	double slope = pt->slope(pt2);
	cout << "Distance: " << dist <<", Alt: " << alt << ", Azim: " << azim << ", Slope: " << slope << "\n" << flush;
	pt->update(10.0, 10.0, 10.0);
	pt2->update(20.0, 20.0, 20.0);
	dist = pt->distance(pt2);
	alt = pt->altitd(pt2);
	azim = pt->azimth(pt2);
	slope = pt->slope(pt2);
	cout << "Distance: " << dist << ", Alt: " << alt << ", Azim: " << azim << ", Slope: " << slope << "\n" << flush;
	pt->print();
	CRS *crs = new CRS();
	crs->transform(pt);
	crs->transform(pt2);
	delete crs;
	pt->print();
	pt2->print();
	dist = pt->distance(pt2);
	alt = pt->altitd(pt2);
	azim = pt->azimth(pt2);
	slope = pt->slope(pt2);
	cout << "Distance: " << dist << ", Alt: " << alt << ", Azim: " << azim << ", Slope: " << slope << "\n" << flush;
	delete pt;
	delete pt2;

	return 0;
}
