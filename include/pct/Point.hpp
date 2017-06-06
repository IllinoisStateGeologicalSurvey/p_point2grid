#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>

using namespace std;

struct point
{
	double x;
	double y;
	double z;
	point();
	point(double x, double y, double z): x(x), y(y), z(z) {};
	point(const point &pt);
	~point();
	void clear();
	//int orientation(point pt);
	void update(double x, double y, double z);
	void print() const;
	//bool operator < (const point *pt2);
	//bool operator > (const point *pt2);
	
	int  orientation(const point *pt);
	double magntd() const;
	double dot(const point *pt) const;
	double distance(const point *pt) const;
	double altitd(const point *pt2) const;
	double azimth(const point *pt2) const;
	double slope(const point *pt2) const;
	void randomize();
};

bool compareX(const point &pt1, const point &pt2);
bool compareY(const point &pt1, const point &pt2);
bool compareZ(const point &pt1, const point &pt2);


//bool operator > (const point *pt1, const point *pt2);


#endif 
