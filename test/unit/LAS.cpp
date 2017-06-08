#include <points2grid/Global.hpp>
#include <points2grid/lasfile.hpp>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <float.h>
#include "Point.hpp"
#include "Grid.hpp"
#include "DTypes.h"

int main(int argc, char **argv) 
{
	int buffer_size = 10000;
	char inputFile[1024] = "/home/ncasler/app-dev/point-cloud-tree/data/1000_908.las";
	char inputFile2[1024] = "/home/ncasler/app-dev/point-cloud-tree/data/998_988.las";
	int epsg_code = 0;
	int point_count = 0;
	double min_x, min_y, min_z, max_x, max_y, max_z;
	DType d_type = DT_Float32;
	int np = 10;
	min_x = DBL_MAX;
	min_y = DBL_MAX;
	min_z = DBL_MAX;
	max_x = -DBL_MAX;
	max_y = -DBL_MAX;
	max_z = -DBL_MAX;
	printf("Input file: %s\n", inputFile);
	las_file las1;
	las_file las2;
	las1.open(inputFile);
	las2.open(inputFile2);
	double resX = 3.0;
	double resY = 3.0;
	min_x = fmin(las1.minimums()[0], las2.minimums()[0]);
	min_y = fmin(las1.minimums()[1], las2.minimums()[1]);
	min_z = fmin(las1.minimums()[2], las2.minimums()[2]);
	max_x = fmax(las1.maximums()[0], las2.maximums()[0]);
	max_y = fmax(las1.maximums()[1], las2.maximums()[1]);
	max_z = fmax(las1.maximums()[2], las2.maximums()[2]);
	printf("MIN:(%f, %f),Max (%f,%f)\n", min_x, min_y, max_x, max_y);
	int cols = ceil((max_x - min_x) / resX);
	int rows = ceil((max_y - min_y) / resY);
	struct point origin(min_x, max_y, min_z);
	struct grid *gridTest = new grid(&origin, cols, rows, d_type, resX, resY);
	int size = gridTest->getSize() / 1000000;
	printf("Grid will have %ix%i dims and will take %i GB\n", cols, rows, size);
	las1.close();
	las2.close();
	// Optimize block size
	double aspect_ratio = (float)cols / (float)rows;
	int cellCount = cols * rows;
	int blockMax = 2000000;
	int b_dim = ceil(sqrt(blockMax));
	int nb = ceil((float)cellCount/(float)blockMax);
	int b_col = ceil((float)cols / (float)b_dim);
	int b_row = ceil((float)rows / (float)b_dim);
	printf("Block_dim: %i, b_col: %i, b_row: %i\n", b_dim, b_col, b_row);
	printf("t_Cell: %i, Proportion: %f, Block Count: %i", cellCount, aspect_ratio, nb);
	delete(gridTest);
	return 0;
}
