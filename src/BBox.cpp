#include <iostream>
#include <vector>
#include "Point.hpp"
#include "BBox.hpp"
/*
BBox::BBox() {
	struct Point min;
	struct Point max;
};

BBox::BBox(struct Point * _min, struct Point *_max) {
	struct Point min(_min->x, _min->y);
	struct Point max(_max->x, _max->y);
};
*/
/* Copy Constructor */
/*BBox::BBox(const BBox& box) {
	struct Point min(box.min);
	struct Point max(box.max);
}*/

BBox::~BBox() {
	min.update(0.0,0.0,0.0);
	max.update(0.0,0.0,0.0);
	//delete min;
	//delete max;
};

void BBox::update(struct Point *_min, struct Point *_max) {
	min.update(_min->x, _min->y, _min->z);
	max.update(_max->x, _max->y, _min->z);
};

void BBox::updateMin(double _x, double _y, double _z) {
	min.update(_x, _y, _z);
};

void BBox::updateMax(double _x, double _y, double _z) {
	max.update(_x, _y, _z);
};

bool BBox::validate() {
	if ((max.x > min.x && max.y > min.y) && max.z > min.z) {
		return true;
	}
	else {
		return false;
	}
}

struct Point* BBox::centroid() {
	//min.print();
	//max.print();
	
	double midX = ((max.x - min.x)/2) + min.x;
	double midY = ((max.y - min.y)/2) + min.y;
	double midZ = ((max.z - min.z)/2) + min.z;
	struct Point *center = new Point(midX, midY, midZ);
	return center;
}

double BBox::area() {
	double vert =  max.y - min.y;
	double width =  max.x - min.x;
	double area = vert * width;
	return area;
}




vector<BBox*> BBox::subdivide() {
	struct Point* center = centroid();
	//cout << "Creating centroid" << endl;
	//center->print();
	vector<BBox*> splits;
	//NW Box
	
	struct Point _min(min.x, center->y, min.z);
	struct Point _max(center->x, max.y, max.z);
	BBox box;
	//cout << "NW" << endl;
	splits.push_back(new BBox(&_min, &_max));
	//cout << "Subdivision Area: " << area() << endl;
	// NE Box
	//cout << "NE" << endl;
	_min.update(center->x, center->y, min.z);
	_max.update(max.x, max.y, max.z);
	splits.push_back(new BBox(&_min, &_max));
	// SW Box
	//cout << "SW" << endl;
	_min.update(min.x, min.y, min.z);
	_max.update(center->x, center->y, max.z);
	splits.push_back(new BBox(&_min, &_max));
	// SE box
	//cout << "SE" << endl;
	_min.update(center->x, min.y, min.z);
	_max.update(max.x, center->y, max.z);
	splits.push_back(new BBox(&_min, &_max));
	// Clean up
	for (vector<BBox*>::iterator it = splits.begin(); it != splits.end(); ++it) {
	//	(*it)->print();
	}
	//cout << "Cleaning up"<< endl;
	//delete box;
	//delete _min;
	//delete _max;
	return splits;
}

void BBox::print() const {
	cout << "   ---------";
	max.print();
	cout << endl;
	int i = 0;
	for (i = 0; i < 5; i++) {
		cout << "   |                              |" << endl;
	}
	cout << "   ";
	min.print();
	cout << "--------" << endl;
};
