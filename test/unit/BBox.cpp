#include <iostream>
#include <vector>
#include "Point.hpp"
#include "BBox.hpp"

int main(int argc, char * argv[]) {
	struct point *min = new point(0.0, 10.0, -30.0);
	struct point *max = new point(20.0, 30.0, 30.0);
	BBox* box = new BBox(min, max);
	box->print();

	min->update(30.0, 40.0, -100.0);
	max->update(60.0, 70.0, 50.0);
	box->update(min, max);

	box->print();
	
	vector<BBox*> children = box->subdivide();
	for (vector<BBox*>::iterator it = children.begin(); it != children.end(); ++it) {
		(*it)->print();
	}
	for (vector<BBox*>::iterator it = children.begin(); it != children.end(); ++it) {
		delete (*it);
	}
	children.clear();
	delete box;
	delete min;
	delete max;

	return 0;
}
