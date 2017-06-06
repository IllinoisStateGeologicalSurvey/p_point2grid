#include <iostream>
#include <iomanip>
#include <vector>
#include "Point.hpp"
#include "util.hpp"
#include <math.h>

using namespace std;

point::point() {
	x = 0;
	y = 0;
	z = 0;
}

point::~point() {
	x = 0;
	y = 0;
	z = 0;
}

/* Copy constructor */
point::point(const struct point &pt) {
	x = pt.x;
	y = pt.y;
	z = pt.z;
};

void point::update(double _x, double _y, double _z) {
	x = _x;
	y = _y;
	z = _z;
}

void point::clear() {
	x = 0.0;
	y = 0.0;
	z = 0.0;
};


void point::print() const {
	// Set desired precision
	cout << fixed;
	cout << setprecision(2);
	cout << "(" << x << ", " << y << ", " << z << ")";
	cout << "\n" << flush;
}

// Currently checks for that relation holds true for all dimensions
// May not be best definition, but currently makes bbox validation 
// easier
//bool point::operator < (const struct point *pt2)
//	{
//		return (x < pt2->x
//				&& y < pt2->y);
//}


//bool point::operator > (const struct point *pt2) {
//	return (x > pt2->x
//			&& y > pt2->y);
//}
// Get vector magnitude(or length from origin)
double point::magntd() const {
	return sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
}
// Return dot product of points
double point::dot(const struct point *pt) const {
	return (x * pt->x + y * pt->y + z * pt->z);
}


double point::distance(const struct point *pt2) const {
	double dx = x - pt2->x;
	double dy = y - pt2->y;
	double dz = z - pt2->z;
	return sqrt(pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0));
	//return (dot(pt2)/(magntd() * pt2->magntd()));
}

double point::altitd(const struct point *pt2) const {
	double dx = x - pt2->x;
	double dy = y - pt2->y;
	double dz = z - pt2->z;
	return to_degrees(atan2(dy, sqrt(dx*x + dz*z)));
}

double point::slope(const struct point *pt2) const {
	double dist = distance(pt2);
	double dz = pt2->z - z;
	return (dz / dist);
}

double point::azimth(const struct point *pt2) const {
	double dx = x - pt2->x;
	double dy = y - pt2->y;
	double dz = z - pt2->z;
	return to_degrees(atan2(-dx, -dz));
}

// Returns a flag to denote the orienation in terms of quadrants
//  Orienation follows nw in a z-order such as:
//
//              -------------
//              |  0  |  1  |
//              |-----------|
//              |  2  |  3  |
//              -------------
int point::orientation(const struct point *pt) {
	int flag = 0;
	if (pt->x > x) {
		flag++;
	}
	if (pt->y > y) {
		flag += 2;
	}
	return flag;
}



void point::randomize() {
	x = randomVal(360, -180);
	y = randomVal(180, -90);
	z = randomVal(12500, -12500);
}


// Methods to compare for sorting
bool compareX(const point &pt1, const point &pt2) {
	return pt1.x < pt2.x;
}

bool compareY(const point &pt1, const point &pt2) {
	return pt1.y < pt2.y;
}

bool compareZ(const point &pt1, const point &pt2) {
	return pt1.z < pt2.z;
}
