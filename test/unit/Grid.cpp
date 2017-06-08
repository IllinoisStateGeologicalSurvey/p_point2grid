#include <iostream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include "pct/Point.hpp"
#include "pct/Grid.hpp"
#include "pct/DTypes.hpp"
#include <stdio.h>
using namespace std;

int main(int argc, char* argv[]) {
	struct Point origin;
	DType d_type = DT_Float32;
	origin.randomize();
	origin.update(-142.0, 42.0, 1000);
	origin.print();
	struct Point test;
	test.randomize();
	test.update(-140.0, 40.4, 1000);
	test.print();
	int cols = 1000;
	int rows = 1000;
	double resX = 1.0;
	double resY = 1.0;
	int idx[2] = {0,0};
	struct Grid *grid1 = new Grid(&origin, cols, rows, d_type, resX, resY);
	cout << "Grid Data type: " << GetDataTypeName(d_type) << ", cell Count: " << cols * rows << ", Total Size: " << grid1->getSize() << "\n";
	
	double Max[2] = {grid1->getMaxX(), origin.y};
	double Min[2] = {origin.x, grid1->getMinY()};
	cout << "Min: " << Min[0] << ", " << Min[1] << "], Max: " << Max[0] << ", " << Max[1] << "\n";
	
	cout << "Allocating cells\n";
	grid1->alloc();
	int i, j = 0;
	// Attempt to set values for grid
	for (i =0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			test.update(j, i, j+i);
			grid1->set(&test);
		}
	}
	Pixel* pix = grid1->get(cols-1, rows-1);
	printf("Pixel at [%i,%i] is %f", cols-1, rows-1, pix->sum);
	char outPath[1024] = "/gpfs_scratch/ncasler/data/tmp/grid.tif";
	int epsg = 3443;
	grid1->write(&outPath[0], epsg);

	int sizeGB = grid1->getSize() / 1000000;
	cout << "Successfully allocated: " << sizeGB << " GB of memory\n";
	grid1->dealloc();
	grid1->getCell(&test, &idx[0]);
	cout << "Idx: " << idx[0] << ", " << idx[1] << "\n";
	
	delete(grid1);
	return 0;

}
