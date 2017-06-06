#ifndef BBOX_HPP
#define BBOX_HPP

#include <iostream>
#include <vector>
#include "Point.hpp"

using namespace std;

class BBox
{
	private:
		struct point min;
		struct point max;
	public:
		BBox(): min(0.0, 0.0, 0.0), max(0.0, 0.0, 0.0){};
		BBox(struct point * _min, struct point * _max): min(*_min), max(*_max){};
		BBox(const BBox& box): min(box.min), max(box.max){};
		~BBox();
		void update(struct point * _min, struct point * _max);
		void updateMin(double _x, double _y, double _z);
		void updateMax(double _x, double _y, double _z);
		struct point* centroid();
		double area();
		vector<BBox*> subdivide();
		bool validate();
		void print() const;
};



#endif
