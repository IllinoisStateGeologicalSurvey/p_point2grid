#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>

using namespace std;

typedef struct Point
{
	double x;
	double y;
	double z;
	Point();
	Point(double x, double y, double z): x(x), y(y), z(z) {};
	Point(const Point &pt);
	~Point();
	void clear();
	//int orientation(point pt);
	void update(double x, double y, double z);
	void print() const;
	//bool operator < (const point *pt2);
	//bool operator > (const point *pt2);
	
	int  orientation(const Point *pt);
	double magntd() const;
	double dot(const Point *pt) const;
	double distance(const Point *pt) const;
	double altitd(const Point *pt2) const;
	double azimth(const Point *pt2) const;
	double slope(const Point *pt2) const;
	void randomize();
} Point;

bool compareX(const Point &pt1, const Point &pt2);
bool compareY(const Point &pt1, const Point &pt2);
bool compareZ(const Point &pt1, const Point &pt2);


//bool operator > (const point *pt1, const point *pt2);


#endif 
