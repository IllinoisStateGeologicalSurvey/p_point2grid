#include <iostream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include "Point.hpp"
#include "BBox.hpp"
#include "PointBucket.hpp"
#include "QuadNode.hpp"
#include "util.hpp"
using namespace std;

int main(int argc, char* argv[]) 
{
	struct point* pt1 = new point(-180, -90, -12500);
	struct point* pt2 = new point(180, 90, 12500);
	//pt1->randomize();
	//pt2->randomize();
	double vm, rss;
	BBox* box = new BBox();
	if (pt1 < pt2) {
		box->update(pt1, pt2);
	} else {
		box->update(pt2, pt1);
	}
	//box->print();
	struct QuadNode* node = new QuadNode(1000000, 0, box, "/tmp");
	node->bounds->print();
	//for (int i = 0; i < 1000; i++) {
	//	pt1->randomize();
	//	node->insert(pt1);
	//}
	//box->print();
	size_t bucket_size = node->bucket->count;
	//cout << "Bucket size: "<< bucket_size << endl;
	// Note this fails somewhere between 647000 and 648000
	//for (int i = 0; i < 647960; i++) {
	//	pt1->randomize();
		//pt1->print();
	//	node->insert(*pt1);
	//}
	for (size_t i = 0; i < 200000000; i++) {
		//cout << "PT: " << i << endl;
		pt1->randomize();
		node->insert(pt1);
		if (i % 10000 == 0) {
			process_mem_usage(vm, rss);
			cout << "VM: " << vm << "; RSS: " << rss << "Count: " << i << endl;
		}
	}
	//pt1->randomize();
	//node->insert(*pt1);
	//node->bucket->print();
	//cout << node->bucket->pts->size() << endl;
	cout << "Size ofL Point: "<< sizeof(*pt1) << " sizeof Node: "<< sizeof(node) << endl;
	delete node;
	delete box;
	delete pt2;
	delete pt1;
};
