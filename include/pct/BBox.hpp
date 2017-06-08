#ifndef BBOX_HPP
#define BBOX_HPP

#include <iostream>
#include <vector>
#include "Point.hpp"

using namespace std;

class BBox
{
	private:
		struct Point min;
		struct Point max;
	public:
		BBox(): min(0.0, 0.0, 0.0), max(0.0, 0.0, 0.0){};
		BBox(struct Point * _min, struct Point * _max): min(*_min), max(*_max){};
		BBox(const BBox& box): min(box.min), max(box.max){};
		~BBox();
		void update(struct Point * _min, struct Point * _max);
		void updateMin(double _x, double _y, double _z);
		void updateMax(double _x, double _y, double _z);
		struct Point* centroid();
		double area();
		vector<BBox*> subdivide();
		bool validate();
		int intersect(const BBox* other);
		void print() const;
};



#endif
