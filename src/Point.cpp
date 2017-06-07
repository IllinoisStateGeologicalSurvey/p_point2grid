#include <iostream>
#include <iomanip>
#include <vector>
#include "pct/Point.hpp"
#include "pct/util.hpp"
#include <math.h>

using namespace std;

Point::Point() {
	x = 0;
	y = 0;
	z = 0;
}

Point::~Point() {
	x = 0;
	y = 0;
	z = 0;
}

/* Copy constructor */
Point::Point(const struct Point &pt) {
	x = pt.x;
	y = pt.y;
	z = pt.z;
};

void Point::update(double _x, double _y, double _z) {
	x = _x;
	y = _y;
	z = _z;
}

void Point::clear() {
	x = 0.0;
	y = 0.0;
	z = 0.0;
};


void Point::print() const {
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
double Point::magntd() const {
	return sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
}
// Return dot product of points
double Point::dot(const struct Point *pt) const {
	return (x * pt->x + y * pt->y + z * pt->z);
}


double Point::distance(const struct Point *pt2) const {
	double dx = x - pt2->x;
	double dy = y - pt2->y;
	double dz = z - pt2->z;
	return sqrt(pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0));
	//return (dot(pt2)/(magntd() * pt2->magntd()));
}

double Point::altitd(const struct Point *pt2) const {
	double dx = x - pt2->x;
	double dy = y - pt2->y;
	double dz = z - pt2->z;
	return to_degrees(atan2(dy, sqrt(dx*x + dz*z)));
}

double Point::slope(const struct Point *pt2) const {
	double dist = distance(pt2);
	double dz = pt2->z - z;
	return (dz / dist);
}

double Point::azimth(const struct Point *pt2) const {
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
int Point::orientation(const struct Point *pt) {
	int flag = 0;
	if (pt->x > x) {
		flag++;
	}
	if (pt->y > y) {
		flag += 2;
	}
	return flag;
}



void Point::randomize() {
	x = randomVal(360, -180);
	y = randomVal(180, -90);
	z = randomVal(12500, -12500);
}


// Methods to compare for sorting
bool compareX(const Point &pt1, const Point &pt2) {
	return pt1.x < pt2.x;
}

bool compareY(const Point &pt1, const Point &pt2) {
	return pt1.y < pt2.y;
}

bool compareZ(const Point &pt1, const Point &pt2) {
	return pt1.z < pt2.z;
}
